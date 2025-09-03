/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FaceCellWave.H"
#include "wallDistAddressing.H"
#include "wallPointAddressing.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "patchDistMethod.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDistAddressing, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::wallDistAddressing::getValues
(
    const FaceCellWave<wallPointAddressing>& wave,
    const List<wallPointAddressing>& allCellInfo,
    const List<wallPointAddressing>& allFaceInfo,
    volScalarField& y
) const
{
    // Extract volField. See patchWave::getValues

    label nIllegal = 0;
    forAll(allCellInfo, celli)
    {
        const scalar dist = allCellInfo[celli].distSqr();

        if (allCellInfo[celli].valid(wave.data()))
        {
            y[celli] = Foam::sqrt(dist);
        }
        else
        {
            y[celli] = dist;
            ++nIllegal;
        }
    }

    // Copy boundary values
    auto& bfld = y.boundaryFieldRef();

    for (auto& pfld : bfld)
    {
        scalarField patchField(pfld.size(), Zero);

        forAll(pfld, patchFacei)
        {
            const label meshFacei = pfld.patch().start() + patchFacei;
            const scalar dist = allFaceInfo[meshFacei].distSqr();

            if (allFaceInfo[meshFacei].valid(wave.data()))
            {
                // Adding SMALL to avoid problems with /0 in the turbulence
                // models
                patchField[patchFacei] = Foam::sqrt(dist) + SMALL;
            }
            else
            {
                patchField[patchFacei] = dist;
                ++nIllegal;
            }
        }
        pfld = patchField;
    }

    return nIllegal;
}


void Foam::wallDistAddressing::addItem
(
    const label item,
    const labelPair& data,
    label& untransformi,
    label& transformi,
    labelPairList& transformedWallInfo
)
{
    const auto& gt = mesh_.globalData().globalTransforms();

    if (gt.transformIndex(data) == gt.nullTransformIndex())
    {
        // Store item and global index of untransformed data
        const label proci = gt.processor(data);
        const label index = gt.index(data);
        untransformedItems_[untransformi] = item;
        untransformedSlots_[untransformi++] =
            globalWallsPtr_().toGlobal(proci, index);
    }
    else
    {
        // Store item and transformed data
        transformedItems_[transformi] = item;
        transformedWallInfo[transformi++] = data;
    }
}


void Foam::wallDistAddressing::correct(volScalarField& y)
{
    y = dimensionedScalar("yWall", dimLength, GREAT);

    const auto& C = mesh_.C();
    const auto& gt = mesh_.globalData().globalTransforms();

    label nWalls = 0;
    for (const label patchi : patchIDs_)
    {
        nWalls += C.boundaryField()[patchi].size();
    }

    globalWallsPtr_.reset(new globalIndex(nWalls));
    const globalIndex& globalWalls = globalWallsPtr_();

    DynamicList<label> seedFaces(nWalls);
    DynamicList<wallPointAddressing> seedInfo(nWalls);


    nWalls = 0;
    for (const label patchi : patchIDs_)
    {
        const auto& fc = C.boundaryField()[patchi];
        const auto areaFraction(fc.patch().patch().areaFraction());

        forAll(fc, patchFacei)
        {
            if
            (
                !areaFraction
             || (areaFraction()[patchFacei] > 0.5)    // mostly wall
            )
            {
                seedFaces.append(fc.patch().start()+patchFacei);
                seedInfo.append
                (
                    wallPointAddressing
                    (
                        fc[patchFacei],
                        gt.encode
                        (
                            Pstream::myProcNo(),
                            nWalls,
                            gt.nullTransformIndex()
                        ),
                        scalar(0.0)
                    )
                );
            }
            nWalls++;
        }
    }

    List<wallPointAddressing> allCellInfo(mesh_.nCells());
    // Initialise the passive data (since no mesh reference in the constructor
    // of wallPointAddressing)
    forAll(allCellInfo, celli)
    {
        allCellInfo[celli].data() = gt.encode
        (
            Pstream::myProcNo(),
            celli,
            gt.nullTransformIndex()
        );
    }
    List<wallPointAddressing> allFaceInfo(mesh_.nFaces());
    forAll(allFaceInfo, facei)
    {
        allFaceInfo[facei].data() = gt.encode
        (
            Pstream::myProcNo(),
            facei,
            gt.nullTransformIndex()
        );
    }

    // Propagate information inwards
    FaceCellWave<wallPointAddressing> wave
    (
        mesh_,
        seedFaces,
        seedInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1
    );


    // Extract wall distance
    // ~~~~~~~~~~~~~~~~~~~~~

    getValues(wave, allCellInfo, allFaceInfo, y);


    const labelHashSet patchSet(patchIDs_);

    // Correct wall cells for true distance
    // - store for corrected cells the nearest wall face
    // - update distance field (y) for these cells
    Map<label> cellToWallFace;
    if (correctWalls_)
    {
        cellToWallFace.reserve(nWalls);

        if (cellDistFuncs::useCombinedWallPatch)
        {
            // Correct across multiple patches
            correctBoundaryCells
            (
                patchIDs_,
                true,           // do point-connected cells as well
                y,
                cellToWallFace
            );
        }
        else
        {
            // Optional backwards compatibility
            correctBoundaryFaceCells
            (
                patchSet,
                y,
                cellToWallFace
            );
            correctBoundaryPointCells
            (
                patchSet,
                y,
                cellToWallFace
            );
        }
    }

    // Make sure boundary values are up to date
    y.correctBoundaryConditions();


    // Extract all addressing
    // ~~~~~~~~~~~~~~~~~~~~~~
    //  - untransformed : -local cell/boundary, -nearest (global) wall
    //  - transformed : -local cell/boundary, -nearest (global) wall

    label untransformi = 0;
    untransformedItems_.resize(mesh_.nCells()+mesh_.nBoundaryFaces());
    untransformedSlots_.resize(untransformedItems_.size());

    label transformi = 0;
    transformedItems_.resize(mesh_.nCells()+mesh_.nBoundaryFaces());
    labelPairList transformedWallInfo(transformedItems_.size());

    for (label celli = 0; celli < mesh_.nCells(); celli++)
    {
        const auto cellFnd = cellToWallFace.find(celli);
        if (cellFnd)
        {
            const label wallFacei = cellFnd();
            const auto& data = allFaceInfo[wallFacei].data();
            addItem(celli, data, untransformi, transformi, transformedWallInfo);
        }
        else if (allCellInfo[celli].valid(wave.data()))
        {
            const auto& data = allCellInfo[celli].data();
            addItem(celli, data, untransformi, transformi, transformedWallInfo);
        }
    }
    untransformedPatchStarts_.resize(mesh_.boundary().size()+1);
    transformedPatchStarts_.resize(mesh_.boundary().size()+1);
    for (const auto& pp : mesh_.boundary())
    {
        untransformedPatchStarts_[pp.index()] = untransformi;
        transformedPatchStarts_[pp.index()] = transformi;

        forAll(pp, patchFacei)
        {
            const label facei = pp.start()+patchFacei;
            if (allFaceInfo[facei].valid(wave.data()))
            {
                const auto& data = allFaceInfo[facei].data();
                addItem
                (
                    facei,
                    data,
                    untransformi,
                    transformi,
                    transformedWallInfo
                );
            }
        }
    }

    untransformedItems_.resize(untransformi);
    untransformedSlots_.resize(untransformi);
    untransformedPatchStarts_.back() = untransformi;
    transformedItems_.resize(transformi);
    transformedWallInfo.resize(transformi);
    transformedPatchStarts_.back() = transformi;

    if (debug)
    {
        Pout<< typeName
            << " : untransformed:" << untransformi
            << " transformed:" << transformi
            << " out of nWalls:" << nWalls
            << " untransformStart:" << flatOutput(untransformedPatchStarts_)
            << " transformStart:" << flatOutput(transformedPatchStarts_)
            << endl;
    }

    List<Map<label>> compactMap;
    mapPtr_.reset
    (
        new mapDistribute
        (
            globalWalls,            // globalIndex
            untransformedSlots_,    // untransformedElements

            gt,                     // globalIndexAndTransform
            transformedWallInfo,    // transformedElements
            transformedSlots_,      // transformedIndices
            compactMap
        )
    );

    if (debug & 2)
    {
        // Use the obtained map to write nearest
        {
            OBJstream os(mesh_.polyMesh::path()/"nearest.obj");
            Pout<< typeName
                << " : writing line from cell/face to nearest wall face to "
                << os.name() << endl;

            // Seed all wall faces with centroid and override all other
            // cells/faces with that of nearest wall
            volVectorField wallCentre(mesh_.C());
            this->map(wallCentre, mapDistribute::transformPosition());

            // Seed all wall faces with normal and override all other
            // cells/faces with that of nearest wall
            volVectorField wallNormal
            (
                IOobject
                (
                    "n" & patchTypeName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobjectOption::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(dimless, Zero),
                patchDistMethod::patchTypes<vector>(mesh_, patchSet)
            );
            // Calculate normal
            for (const label patchi : patchIDs_)
            {
                auto& pnf = wallNormal.boundaryFieldRef()[patchi];
                pnf == pnf.patch().nf();
            }
            this->map(wallNormal, mapDistribute::transform());

            forAll(wallCentre, celli)
            {
                const point& cc = mesh_.C()[celli];
                const point& wallC = wallCentre[celli];
                const vector& wallN = wallNormal[celli];
                os.write(linePointRef(cc, wallC), wallN, wallN);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDistAddressing::wallDistAddressing
(
    const fvMesh& mesh,
    const bool correctWalls,
    const label updateInterval
)
:
    // Register as "wallDistAddressing"
    MeshObject_type(mesh),
    cellDistFuncs(mesh),
    patchIDs_(mesh.boundaryMesh().findPatchIDs<wallPolyPatch>().sortedToc()),
    patchTypeName_("wall"),
    updateInterval_(updateInterval),
    correctWalls_(correctWalls),
    requireUpdate_(true),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)    // Same as wallDistData
    )
{
    movePoints();
}


Foam::wallDistAddressing::wallDistAddressing
(
    const word& patchTypeName,
    const fvMesh& mesh,
    const labelList& patchIDs,
    const bool correctWalls,
    const label updateInterval
)
:
    MeshObject_type(patchTypeName, mesh),
    cellDistFuncs(mesh),
    patchIDs_(patchIDs),
    patchTypeName_(patchTypeName),
    updateInterval_(updateInterval),
    correctWalls_(correctWalls),
    requireUpdate_(true),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)    // Same as wallDistData
    )
{
    movePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDistAddressing::~wallDistAddressing()
{}  // mapDistribute was forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallDistAddressing::movePoints()
{
    if
    (
        (updateInterval_ > 0)
     && ((mesh_.time().timeIndex() % updateInterval_) == 0)
    )
    {
        requireUpdate_ = true;
    }

    if (requireUpdate_)
    {
        DebugInfo<< "Updating wall distance" << endl;

        requireUpdate_ = false;

        correct(y_);

        return true;
    }

    return false;
}


void Foam::wallDistAddressing::updateMesh(const mapPolyMesh& mpm)
{
    // Force update if performing topology change
    // Note: needed?
    // - field would have been mapped, so if using updateInterval option (!= 1)
    //   live with error associated of not updating and use mapped values?
    requireUpdate_ = true;
    movePoints();
}


// ************************************************************************* //
