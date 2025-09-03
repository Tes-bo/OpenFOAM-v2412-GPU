/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "parallelFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "primitiveMeshTools.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parallelFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        parallelFvGeometryScheme,
        dict
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parallelFvGeometryScheme::adjustGeometry
(
    pointField& faceCentres,
    vectorField& faceAreas
) const
{
    pointField syncedBCentres;
    syncedBCentres = SubField<vector>
    (
        faceCentres,
        mesh_.nBoundaryFaces(),
        mesh_.nInternalFaces()
    );
    syncTools::swapBoundaryFaceList
    (
        mesh_,
        syncedBCentres
    );

    vectorField syncedBAreas;
    syncedBAreas = SubField<vector>
    (
        faceAreas,
        mesh_.nBoundaryFaces(),
        mesh_.nInternalFaces()
    );
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        syncedBAreas,
        eqOp<vector>(),
        transformOriented()
    );

    const auto& pbm = mesh_.boundaryMesh();
    for (const auto& pp : pbm)
    {
        const auto* ppp = isA<coupledPolyPatch>(pp);

        //if (ppp)
        //{
        //    Pout<< "For patch:" << ppp->name()
        //        << " size:" << ppp->size()
        //        << endl;
        //    forAll(*ppp, i)
        //    {
        //        const label facei = ppp->start()+i;
        //        Pout<< "    Face:" << facei << nl
        //            << "    meshFc:" << faceCentres[facei] << nl
        //            << "    meshFa:" << faceAreas[facei] << nl
        //            << endl;
        //    }
        //}

        if (ppp && !ppp->owner())
        {
            SubField<point> patchFc
            (
                faceCentres,
                ppp->size(),
                ppp->start()
            );
            patchFc = SubField<vector>
            (
                syncedBCentres,
                ppp->size(),
                ppp->offset()
            );

            SubField<vector> patchArea
            (
                faceAreas,
                ppp->size(),
                ppp->start()
            );
            patchArea = SubField<vector>
            (
                syncedBAreas,
                ppp->size(),
                ppp->offset()
            );
        }
    }
}


void Foam::parallelFvGeometryScheme::adjustGeometry()
{
    // Swap face centres and areas
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField faceCentres(mesh_.faceCentres());
    vectorField faceAreas(mesh_.faceAreas());
    {
        pointField syncedBCentres;
        syncedBCentres = SubField<vector>
        (
            mesh_.faceCentres(),
            mesh_.nBoundaryFaces(),
            mesh_.nInternalFaces()
        );
        syncTools::swapBoundaryFaceList
        (
            mesh_,
            syncedBCentres
        );

        vectorField syncedBAreas;
        syncedBAreas = SubField<vector>
        (
            mesh_.faceAreas(),
            mesh_.nBoundaryFaces(),
            mesh_.nInternalFaces()
        );
        syncTools::syncBoundaryFaceList
        (
            mesh_,
            syncedBAreas,
            eqOp<vector>(),
            transformOriented()
        );

        const auto& pbm = mesh_.boundaryMesh();
        for (const auto& pp : pbm)
        {
            const auto* ppp = isA<coupledPolyPatch>(pp);

            //if (ppp)
            //{
            //    Pout<< "For patch:" << ppp->name()
            //        << " size:" << ppp->size()
            //        << endl;
            //    forAll(*ppp, i)
            //    {
            //        const label facei = ppp->start()+i;
            //        Pout<< "    Face:" << facei << nl
            //            << "    meshFc:" << faceCentres[facei] << nl
            //            << "    meshFa:" << faceAreas[facei] << nl
            //            << endl;
            //    }
            //}

            if (ppp && !ppp->owner())
            {
                // Delete old cached geometry
                const_cast<coupledPolyPatch&>
                (
                    *ppp
                ).primitivePatch::clearGeom();

                SubField<point> patchFc
                (
                    faceCentres,
                    ppp->size(),
                    ppp->start()
                );
                patchFc = SubField<vector>
                (
                    syncedBCentres,
                    ppp->size(),
                    ppp->offset()
                );

                SubField<vector> patchArea
                (
                    faceAreas,
                    ppp->size(),
                    ppp->start()
                );
                patchArea = SubField<vector>
                (
                    syncedBAreas,
                    ppp->size(),
                    ppp->offset()
                );
            }
        }
    }


    // Force recalculating cell properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField cellCentres(mesh_.nCells());
    scalarField cellVolumes(mesh_.nCells());
    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh_,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes
    );


    label nFaces = 0;
    forAll(faceCentres, facei)
    {
        const auto& oldFc = mesh_.faceCentres()[facei];
        const auto& newFc = faceCentres[facei];
        const auto& oldFa = mesh_.faceAreas()[facei];
        const auto& newFa = faceAreas[facei];

        if (oldFc != newFc || oldFa != newFa)
        {
            nFaces++;

            if (mesh_.isInternalFace(facei))
            {
                WarningInFunction
                    << "Different geometry for internal face:" << facei
                    << "    oldFc:" << oldFc << nl
                    << "    newFc:" << newFc << nl
                    << "    oldFa:" << oldFa << nl
                    << "    newFa:" << newFa << nl
                    << endl;
            }
        }
    }

    label nCells = 0;
    forAll(cellCentres, celli)
    {
        const auto& oldCc = mesh_.cellCentres()[celli];
        const auto& newCc = cellCentres[celli];
        const auto& oldCv = mesh_.cellVolumes()[celli];
        const auto& newCv = cellVolumes[celli];

        if (oldCc != newCc || oldCv != newCv)
        {
            nCells++;
            //Pout<< "cell:" << celli << nl
            //    << "    oldCc:" << oldCc << nl
            //    << "    newCc:" << newCc << nl
            //    << "    oldCv:" << oldCv << nl
            //    << "    newCv:" << newCv << nl
            //    << endl;
        }
    }

    Info<< "parallelFvGeometryScheme::movePoints() : "
        << "adjusted geometry of faces:"
        << returnReduce(nFaces, sumOp<label>())
        << " of cells:"
        << returnReduce(nCells, sumOp<label>())
        << endl;

    const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
    (
        std::move(faceCentres),
        std::move(faceAreas),
        std::move(cellCentres),
        std::move(cellVolumes)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parallelFvGeometryScheme::parallelFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvGeometryScheme(mesh, dict),
    dict_(dict.subOrEmptyDict("geometry"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvGeometryScheme& Foam::parallelFvGeometryScheme::geometry() const
{
    if (!geometryPtr_)
    {
        if (debug)
        {
            Pout<< "parallelFvGeometryScheme::geometry() : "
                << "constructing underlying scheme from " << dict_
                << endl;
        }
        geometryPtr_ = fvGeometryScheme::New
        (
            mesh_,
            dict_,
            basicFvGeometryScheme::typeName
        );
    }
    return geometryPtr_();
}


void Foam::parallelFvGeometryScheme::movePoints()
{
    if (debug)
    {
        Pout<< "parallelFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    // Do any primitive geometry calculation
    const_cast<fvGeometryScheme&>(geometry()).movePoints();

    // Do in-place overriding processor boundary geometry
    adjustGeometry();
}


void Foam::parallelFvGeometryScheme::updateMesh(const mapPolyMesh& mpm)
{
    return const_cast<fvGeometryScheme&>(geometry()).updateMesh(mpm);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::parallelFvGeometryScheme::weights() const
{
    return geometry().weights();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::parallelFvGeometryScheme::deltaCoeffs() const
{
    return geometry().deltaCoeffs();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::parallelFvGeometryScheme::nonOrthDeltaCoeffs() const
{
    return geometry().nonOrthDeltaCoeffs();
}


Foam::tmp<Foam::surfaceVectorField>
Foam::parallelFvGeometryScheme::nonOrthCorrectionVectors() const
{
    return geometry().nonOrthCorrectionVectors();
}


bool Foam::parallelFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    primitiveMeshTools::makeFaceCentresAndAreas
    (
        mesh_,
        points,
        faceCentres,
        cellCentres
    );

    // Parallel consistency
    adjustGeometry(faceCentres, faceAreas);

    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh_,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes
    );

    return true;
}


// ************************************************************************* //
