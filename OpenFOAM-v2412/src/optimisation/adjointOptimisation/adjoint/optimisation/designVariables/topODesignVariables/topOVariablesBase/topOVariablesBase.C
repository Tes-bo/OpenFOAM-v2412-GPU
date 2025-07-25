/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "topOVariablesBase.H"
#include "pointMesh.H"
#include "topOZones.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "wallFvPatch.H"
#include "coupledFvPatch.H"
#include "emptyFvPatch.H"
#include "MeshedSurfaceProxy.H"
#include "surfaceWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topOVariablesBase, 0);
}


// * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * //

Foam::DynamicList<Foam::label> Foam::topOVariablesBase::faceFaces
(
    const label facei
) const
{
    const fvMesh& mesh = zones_.mesh();
    const labelListList& edgeFaces = mesh.edgeFaces();
    DynamicList<label> neighs;
    if (!mesh.isInternalFace(facei))
    {
        const labelList& faceEdges = mesh.faceEdges()[facei];
        for (const label edgei: faceEdges)
        {
            const labelList& edgeIFaces = edgeFaces[edgei];
            for (const label neiFacei : edgeIFaces)
            {
                if (neiFacei != facei && !mesh.isInternalFace(neiFacei))
                {
                    const label patchi =
                        mesh.boundaryMesh().whichPatch(neiFacei);
                    if (!isA<emptyFvPatch>(mesh.boundary()[patchi]))
                    {
                        neighs.push_back(neiFacei);
                    }
                }
            }
        }
    }
    return neighs;
}


bool Foam::topOVariablesBase::addCutBoundaryFaceToIsoline
(
    const label facei,
    const cutFaceIso& cutFace,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs,
    List<DynamicList<label>>& cuttingFacesPerMeshFace
) const
{
    const fvMesh& mesh = zones_.mesh();
    // If this is a boundary face being cut, append it to the iso-surface
    if (!mesh.isInternalFace(facei))
    {
        const label cutPatchi = mesh.boundaryMesh().whichPatch(facei);
        const fvPatch& cutPatch = mesh.boundary()[cutPatchi];
        if (!isA<coupledFvPatch>(cutPatch) && !isA<emptyFvPatch>(cutPatch))
        {
            if
            (
                addCuttingFaceToIsoline
                (
                    cutFace.subFacePoints(),
                    cutPatchi,
                    faceFaces(facei),
                    cuttingFacesPerMeshFace,
                    isoSurfPts,
                    isoSurfFaces,
                    zoneIDs
                )
            )
            {
                cuttingFacesPerMeshFace[facei].push_back
                (
                    isoSurfFaces.size() - 1
                );
                return true;
            }

        }
    }
    return false;
}


bool Foam::topOVariablesBase::addCuttingFaceToIsoline
(
    const DynamicList<point>& facePoints,
    const label nSerialPatches,
    const DynamicList<label>& cellCutFaces,
    const List<DynamicList<label>>& cuttingFacesPerMeshFace,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs
) const
{
    if (facePoints.size() > 2)
    {
        // Check whether any of points of the new iso-surface face are already
        // present in the surface. To reduce the number of comparisons, only
        // points on iso-surface faces cutting the current cell are checked
        labelList uniquePointIDs(facePoints.size(), -1);
        DynamicList<point> uniqueFacePoints(facePoints.size());
        DynamicList<label> uniqueFacePointEdges(facePoints.size());
        label addedPoints = Zero;
        forAll(facePoints, pi)
        {
            bool foundInNei = false;
            for (const label cutMeshFacei : cellCutFaces)
            {
                if
                (
                    isDuplicatePoint
                    (
                        pi,
                        facePoints[pi],
                        cuttingFacesPerMeshFace[cutMeshFacei],
                        isoSurfPts,
                        isoSurfFaces,
                        uniquePointIDs
                    )
                )
                {
                    foundInNei = true;
                    break;
                }
            }
            if (!foundInNei)
            {
                uniquePointIDs[pi] = isoSurfPts.size() + addedPoints;
                ++addedPoints;
                uniqueFacePoints.push_back(facePoints[pi]);
            }
        }

        face isoFace(uniquePointIDs);
        isoSurfPts.append(uniqueFacePoints);
        isoSurfFaces.append(isoFace);
        zoneIDs.append(nSerialPatches);

        return true;
    }
    return false;
}


bool Foam::topOVariablesBase::isDuplicatePoint
(
    const label pointID,
    const vector& pointi,
    const DynamicList<label>& cuttingFaces,
    const DynamicList<point>& isoSurfPts,
    const DynamicList<face>& isoSurfFaces,
    labelList& uniquePointIDs
) const
{
    for (const label cuttingFacei : cuttingFaces)
    {
        const face& cuttingFace = isoSurfFaces[cuttingFacei];
        for (const label neiPi : cuttingFace)
        {
            if (mag(pointi - isoSurfPts[neiPi]) < SMALL)
            {
                uniquePointIDs[pointID] = neiPi;
                return true;
            }
        }
    }
    return false;
}


void Foam::topOVariablesBase::addBoundaryFacesToIsoline
(
    const pointScalarField& pointY,
    const Map<label>& addedFaces,
    const scalar isoValue,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs,
    label& nChangedFaces,
    labelList& changedFaces,
    List<wallPointData<label>>& changedFacesInfo,
    labelList& changedFaceToCutFace,
    List<DynamicList<label>>& cuttingFacesPerMeshFace
)
{
    const fvMesh& mesh = zones_.mesh();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        bool isWall = isA<wallFvPatch>(patch);
        if (!isA<emptyFvPatch>(patch) && !isA<coupledFvPatch>(patch))
        {
            const label start = patch.start();
            forAll(patch, facei)
            {
                const label gFacei = start + facei;
                // On rare occasions, a face value might be negative but all
                // point values are positive. In such a case, the face won't be
                // cut and is seen as a fluid face. Hence, the point values are
                // used to determine the status of the boundary face.
                bool isFluid(true);
                for (const label pti : mesh.faces()[gFacei])
                {
                    isFluid = isFluid && (pointY[pti] >= isoValue);
                }
                if (isFluid && !addedFaces.found(gFacei))
                {
                    // Insert wall faces if they belong to the outer boundary
                    if (isWall)
                    {
                        // Mesh face to changedFace addressing
                        meshFaceToChangedFace_.insert(gFacei, nChangedFaces);
                        // Set origin face for meshWave
                        changedFacesInfo[nChangedFaces] =
                            wallPointData<label>
                            (
                                patch.Cf()[facei],
                                nChangedFaces,
                                0
                            );
                        changedFaces[nChangedFaces] = gFacei;
                        // Origin face-to-cut face
                        changedFaceToCutFace.push_back(isoSurfFaces.size());
                        ++nChangedFaces;
                    }

                    if
                    (
                        addCuttingFaceToIsoline
                        (
                            faces[gFacei].points(points),
                            patchi,
                            faceFaces(gFacei),
                            cuttingFacesPerMeshFace,
                            isoSurfPts,
                            isoSurfFaces,
                            zoneIDs
                        )
                    )
                    {
                        cuttingFacesPerMeshFace[gFacei].push_back
                        (
                            isoSurfFaces.size() - 1
                        );
                    }
                }
            }
        }
    }
}


void Foam::topOVariablesBase::writeSurfaceFiles
(
    const pointField& pts,
    const faceList& faces,
    const labelList& zoneIds,
    const label nSerialPatches
) const
{
    const fvMesh& mesh = zones_.mesh();

    // Write vtp file
    const word timeName = mesh.time().timeName();
    fileName localName("topOIsoSurface" + timeName);
    fileName fname(isoSurfFolder_/localName);
    vtk::surfaceWriter vtkFile(pts, faces, fname);
    vtkFile.beginFile();
    vtkFile.writeGeometry();
    vtkFile.beginCellData(1);
    vtkFile.writeCellData("zoneIds", zoneIds);
    vtkFile.close();

    // It appears that there is no interface to write a multi-zone stl file (?)
    // Follow a process similar to proxySurfaceWriter, but also use the zoneIds.

    // Dummy faceIds.
    // faceMap is built after merging the geometry from all processors, to
    // be based on the global addressing
    labelList faceIds;

    autoPtr<meshedSurf> surf(nullptr);
    if (Pstream::parRun())
    {
        surf.reset
        (
            autoPtr<meshedSurf>::NewFrom<mergedSurf>
            (
                pts, faces, zoneIds, faceIds, surfaceWriter::defaultMergeDim
            )
       );
    }
    else
    {
        surf.reset
        (
            autoPtr<meshedSurf>::NewFrom<meshedSurfRef>
            (
                pts, faces, zoneIds, faceIds
            )
        );
    }

    if (Pstream::master())
    {
        const faceList& surfFaces = surf().faces();
        const labelList& surfZoneIds = surf().zoneIds();

        // Size per zone
        labelList zoneSizes(nSerialPatches + 1, 0);
        forAll(surfFaces, fi)
        {
            ++zoneSizes[surfZoneIds[fi]];
        }

        // Faces passed in previous zones
        labelList cumulZoneSizes(nSerialPatches + 1, 0);
        for (label zi = 1; zi < cumulZoneSizes.size(); ++zi)
        {
            cumulZoneSizes[zi] = cumulZoneSizes[zi - 1] + zoneSizes[zi - 1];
        }

        // Construction of the faceMap
        // ---------------------------
        labelList faceMap(surfFaces.size(), -1);
        // Faces visited so far per zone
        labelList passedFacesPerZone(surfZoneIds.size(), 0);

        forAll(surfFaces, fi)
        {
            const label zi = surfZoneIds[fi];
            faceMap[cumulZoneSizes[zi] + passedFacesPerZone[zi]] = fi;
            ++passedFacesPerZone[zi];
        }

        // Construction of the surfZoneList
        // --------------------------------
        List<surfZone> surfZones(nSerialPatches + 1);
        label zi = 0;
        forAll(mesh.boundary(), pi)
        {
            const fvPatch& patch = mesh.boundary()[pi];
            const word& name = patch.name();
            if (!isA<coupledFvPatch>(patch) && !isA<emptyFvPatch>(patch))
            {
                surfZones[zi] =
                    surfZone(name, zoneSizes[zi], cumulZoneSizes[zi], zi);
                zi++;
            }
        }
        surfZones[zi] =
            surfZone
            (
                "topOPatch",
                zoneSizes[nSerialPatches],
                cumulZoneSizes[nSerialPatches],
                zi
            );
        surfZones.setSize(zi + 1);

        MeshedSurfaceProxy<face>
        (
            surf().points(),
            surfFaces,
            surfZones,
            faceMap
        ).write
        (
            fname + ".stl"
            // BINARY seems to loose the patch names
            //IOstreamOption
            //(
            //    IOstreamOption::streamFormat::ASCII,
            //    IOstreamOption::compressionType::COMPRESSED
            //)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topOVariablesBase::topOVariablesBase
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    localIOdictionary
    (
        IOobject
        (
            "topOVars",
            mesh.time().timeName(),
            fileName("uniform"),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    zones_(mesh, dict),
    isoSurfFolder_
        (mesh.time().globalPath()/"optimisation"/"topOIsoSurfaces"),
    meshFaceToChangedFace_(),
    //changedFacesPerCuttingFace_(),
    surfPoints_(),
    surfFaces_()
{
    mkDir(isoSurfFolder_);
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::topOVariablesBase::sourceTerm
(
    DimensionedField<scalar, volMesh>& field,
    const topOInterpolationFunction& interpolationFunc,
    const scalar betaMax,
    const word& interpolationFieldName
) const
{
    // Interpolate based on beta values
    const scalarField& beta = this->beta().primitiveField();
    interpolationFunc.interpolate(beta, field.field());
    field *= betaMax;
}


void Foam::topOVariablesBase::sourceTermSensitivities
(
    scalarField& sens,
    const topOInterpolationFunction& interpolationFunc,
    const scalar betaMax,
    const word& designVariablesName,
    const word& interpolationFieldName
) const
{
    if (designVariablesName == "beta")
    {
        // Multiply with derivative of the interpolation function
        const scalarField& beta = this->beta().primitiveField();
        sens *= betaMax*interpolationFunc.derivative(beta);
    }
}


void Foam::topOVariablesBase::writeFluidSolidInterface
(
    const volScalarField& indicator,
    const scalar isoValue,
    labelList& changedFaces,
    List<wallPointData<label>>& changedFacesInfo
)
{
    const fvMesh& mesh = zones_.mesh();

    // Current number of mesh faces being cut
    label nChangedFaces(0);

    // Indicator at the mesh points, used to compute the iso-surface
    pointScalarField pointY(volPointInterpolation(mesh).interpolate(indicator));
    if (debug && mesh.time().writeTime())
    {
        pointY.rename("pointY");
        pointY.write();
    }

    // Storage for the cutting face/cell mechanism
    cutCellIso cutCell(mesh, pointY);
    cutFaceIso cutFace(mesh, pointY);
    DynamicList<face> isoSurfFaces(mesh.nFaces()/100);
    DynamicList<point> isoSurfPts(mesh.nPoints()/100);
    DynamicList<label> zoneIDs(mesh.nFaces()/100);

    // Number of patches before processor ones
    label nSerialPatches = mesh.boundaryMesh().nNonProcessor();

    // Map between iso-surface and mesh faces (internal and boundary)
    meshFaceToChangedFace_.clearStorage();

    // Map between changedFace and cutFace
    // Index is the changedFaceID, output is the cutFaceID
    DynamicList<label> changedFaceToCutFace;

    // Per mesh face, its cutting faces
    List<DynamicList<label>> cuttingFacesPerMeshFace(mesh.nFaces());

    // Loop over cells and check whether they are cut by zero iso-surface
    const cellList& cells = mesh.cells();

    forAll(cells, celli)
    {
        if (!cutCell.calcSubCell(celli, isoValue))
        {
            const vector& cuttingCf = cutCell.faceCentre();
            vector cuttingNf = cutCell.faceArea();
            cuttingNf.normalise();
            DynamicList<label> cellCutFaces(cells[celli].size());

            // Loop over faces of the cut cell and check whether they are cut
            // too, setting their distance from the cutting face
            forAll(cells[celli], fi)
            {
                const label facei = cells[celli][fi];
                if (!cutFace.calcSubFace(facei, isoValue))
                {
                    const vector& Cf = mesh.Cf()[facei];
                    const scalar distSqr =
                        magSqr((Cf - cuttingCf) & cuttingNf);
                    const vector lsPoint =
                        Cf + ((cuttingCf - Cf) & cuttingNf)*cuttingNf;
                    cellCutFaces.push_back(facei);

                    if (!meshFaceToChangedFace_.found(facei))
                    {
                        // If the face being cut is a boundary one, part of
                        // it belongs to the iso-surface
                        bool addedToIsoSurf = addCutBoundaryFaceToIsoline
                        (
                            facei,
                            cutFace,
                            isoSurfPts,
                            isoSurfFaces,
                            zoneIDs,
                            cuttingFacesPerMeshFace
                        );
                        if (mesh.isInternalFace(facei) || addedToIsoSurf)
                        {
                            // Mesh face to changedFace addressing
                            meshFaceToChangedFace_.insert(facei, nChangedFaces);
                            // Set origin face for meshWave
                            changedFacesInfo[nChangedFaces] =
                                wallPointData<label>
                                (
                                    lsPoint,
                                    nChangedFaces,
                                    distSqr
                                );
                            changedFaces[nChangedFaces] = facei;
                            // Origin face-to-cut face
                            changedFaceToCutFace.push_back(isoSurfFaces.size());
                            ++nChangedFaces;
                        }
                    }
                    else
                    {
                        label visitedFace = meshFaceToChangedFace_.at(facei);
                        if (distSqr < changedFacesInfo[visitedFace].distSqr())
                        {
                            // Set origin face for meshWave
                            changedFacesInfo[visitedFace] =
                                wallPointData<label>
                                (
                                    lsPoint,
                                    visitedFace,
                                    distSqr
                                );
                            // Origin face-to-cut face
                            changedFaceToCutFace[visitedFace] =
                                isoSurfFaces.size();
                        }
                    }
                }
            }

            if
            (
                addCuttingFaceToIsoline
                (
                    cutCell.facePoints(),
                    nSerialPatches,
                    cellCutFaces,
                    cuttingFacesPerMeshFace,
                    isoSurfPts,
                    isoSurfFaces,
                    zoneIDs
                )
            )
            {
                for (const label facei : cellCutFaces)
                {
                    cuttingFacesPerMeshFace[facei].push_back
                    (
                        isoSurfFaces.size() - 1
                    );
                }
            }
        }
    }

    // Add boundary faces on the initial domain with a positive values on all
    // points to the zero iso-surface
    addBoundaryFacesToIsoline
    (
        pointY,
        meshFaceToChangedFace_,
        isoValue,
        isoSurfPts,
        isoSurfFaces,
        zoneIDs,
        nChangedFaces,
        changedFaces,
        changedFacesInfo,
        changedFaceToCutFace,
        cuttingFacesPerMeshFace
    );

    changedFaces.setSize(nChangedFaces);
    changedFacesInfo.setSize(nChangedFaces);

    // vtp and multi-region stl files holding the current geometry
    surfPoints_.transfer(isoSurfPts);
    surfFaces_.transfer(isoSurfFaces);

    const labelList zoneIds(std::move(zoneIDs));
    writeSurfaceFiles(surfPoints_, surfFaces_, zoneIds, nSerialPatches);

    // Invert changedFace-to-cuttingFace map for the sensitivity computations
    //changedFacesPerCuttingFace_ =
    //    invertOneToMany(surfFaces_.size(), changedFaceToCutFace);

    // Transform origin cut faces to a global numbering

    const label changedFacesOffset =
        globalIndex::calcOffset(changedFaces.size());

    for (auto& info : changedFacesInfo)
    {
        info.data() += changedFacesOffset;
    }
}


// ************************************************************************* //
