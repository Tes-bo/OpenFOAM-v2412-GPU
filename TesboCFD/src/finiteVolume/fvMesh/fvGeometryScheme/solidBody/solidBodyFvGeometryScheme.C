/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "solidBodyFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "primitiveMeshTools.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        solidBodyFvGeometryScheme,
        dict
    );
}


bool Foam::solidBodyFvGeometryScheme::markChanges
(
    const pointField& oldPoints,
    const pointField& currPoints,
    bitSet& isChangedPoint,
    bitSet& isChangedFace,
    bitSet& isChangedCell
) const
{
    isChangedPoint.setSize(oldPoints.size());

    // Check for non-identical points
    forAll(isChangedPoint, pointi)
    {
        isChangedPoint.set(pointi, oldPoints[pointi] != currPoints[pointi]);
    }

    DebugInfo
        << "SBM --- Changed points:"
        << returnReduce(isChangedPoint.count(), sumOp<label>())
        << endl;

    // Quick return if no points have moved
    if (returnReduceAnd(isChangedPoint.none()))
    {
        return false;
    }

    isChangedFace.setSize(mesh_.nFaces());
    isChangedFace = false;

    isChangedCell.setSize(mesh_.nCells());
    isChangedCell = false;

    const auto& pointFaces = mesh_.pointFaces();
    const auto& own = mesh_.faceOwner();
    const auto& nbr = mesh_.faceNeighbour();

    // Identify faces and cells attached to moving points
    for (const label pointi : isChangedPoint)
    {
        for (const auto facei : pointFaces[pointi])
        {
            isChangedFace.set(facei);

            isChangedCell.set(own[facei]);
            if (facei < mesh_.nInternalFaces())
            {
                isChangedCell.set(nbr[facei]);
            }
        }
    }

    DebugInfo
        << "SBM --- Changed cells:"
        << returnReduce(isChangedCell.count(), sumOp<label>())
        << endl;

    return true;
}


void Foam::solidBodyFvGeometryScheme::setMeshMotionData()
{
    if (!cacheInitialised_ || !cacheMotion_)
    {
        DebugInFunction << "Creating cache" << endl;

        changedFaceIDs_.clear();    // used for face areas, meshPhi
        changedPatchIDs_.clear();   // used for meshPhi
        changedCellIDs_.clear();    // used for cell volumes

        const pointField& oldPoints = mesh_.oldPoints();
        const pointField& currPoints = mesh_.points();

        if (oldPoints.size() != currPoints.size())
        {
            FatalErrorInFunction
                << "Old and current points sizes must be the same. "
                << "Old points:" << oldPoints.size()
                << " Current points:" << currPoints.size()
                << abort(FatalError);
        }

        //bitSet changedPoints(oldPoints.size());
        //
        //// Check for non-identical points
        //forAll(changedPoints, pointi)
        //{
        //    changedPoints.set
        //    (
        //        pointi,
        //        oldPoints[pointi] != currPoints[pointi]
        //    );
        //}
        //
        //DebugInfo
        //    << "SBM --- Changed points:"
        //    << returnReduce(changedPoints.count(), sumOp<label>())
        //    << endl;
        //
        //// Quick return if no points have moved
        //if (returnReduceAnd(changedPoints.none()))
        //{
        //    return;
        //}
        //
        //bitSet cellIDs(mesh_.nCells());
        //bitSet faceIDs(mesh_.nFaces());
        //
        //const auto& pointFaces = mesh_.pointFaces();
        //const auto& own = mesh_.faceOwner();
        //const auto& nbr = mesh_.faceNeighbour();
        //
        //// Identify faces and cells attached to moving points
        //for (const label pointi : changedPoints)
        //{
        //    for (const auto facei : pointFaces[pointi])
        //    {
        //        faceIDs.set(facei);
        //
        //        cellIDs.set(own[facei]);
        //        if (facei < mesh_.nInternalFaces())
        //        {
        //            cellIDs.set(nbr[facei]);
        //        }
        //    }
        //}
        //
        //changedCellIDs_ = cellIDs.toc();
        //
        //DebugInfo
        //    << "SBM --- Changed cells:"
        //    << returnReduce(changedCellIDs_.size(), sumOp<label>())
        //    << endl;

        bitSet isChangedPoint;
        bitSet isChangedFace;
        bitSet isChangedCell;
        const bool changed = markChanges
        (
            oldPoints,
            currPoints,
            isChangedPoint,
            isChangedFace,
            isChangedCell
        );

        // Quick return if no points have moved
        if (!changed)
        {
            return;
        }

        changedCellIDs_ = isChangedCell.toc();


        // Construct face and patch ID info

        DynamicList<label> changedFaceIDs(isChangedFace.count());
        DynamicList<label> changedPatchIDs(changedFaceIDs.capacity());
        for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
        {
            if (isChangedFace[facei])
            {
                changedFaceIDs.append(facei);
                changedPatchIDs.append(-1);
            }
        }

        const auto& pbm = mesh_.boundaryMesh();
        for (label patchi = 0; patchi < pbm.size(); ++patchi)
        {
            const polyPatch& pp = pbm[patchi];

            for (const label meshFacei : pp.range())
            {
                if (isChangedFace[meshFacei])
                {
                    changedFaceIDs.append(meshFacei);
                    changedPatchIDs.append(patchi);
                }
            }
        }

        changedFaceIDs_.transfer(changedFaceIDs);
        changedPatchIDs_.transfer(changedPatchIDs);
    }

    cacheInitialised_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyFvGeometryScheme::solidBodyFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicFvGeometryScheme(mesh, dict),
    partialUpdate_(dict.getOrDefault<bool>("partialUpdate", true)),
    cacheMotion_(dict.getOrDefault<bool>("cacheMotion", true)),
    cacheInitialised_(false),
    changedFaceIDs_(),
    changedPatchIDs_(),
    changedCellIDs_()
{
    DebugInFunction
        << "partialUpdate:" << partialUpdate_
        << " cacheMotion:" << cacheMotion_
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidBodyFvGeometryScheme::movePoints()
{
    // Note: not calling fvGeometryScheme::movePoints since we want to perform
    // our own geometry manipulations

    bool haveGeometry =
        mesh_.hasCellCentres()
     && mesh_.hasFaceCentres()
     && mesh_.hasCellVolumes()
     && mesh_.hasFaceAreas();

    if (!haveGeometry)
    {
        DebugInFunction
            << "Creating initial geometry using primitiveMesh::updateGeom"
            << endl;

        const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
        return;
    }

    if (mesh_.moving())
    {
        setMeshMotionData();

        DebugInFunction << "Performing partial meshPhi construction" << endl;

        const pointField& oldPoints = mesh_.oldPoints();
        const pointField& currPoints = mesh_.points();

        if (oldPoints.size() != currPoints.size())
        {
            FatalErrorInFunction
                << "Old and current points sizes must be the same. "
                << "Old points:" << oldPoints.size()
                << " Current points:" << currPoints.size()
                << abort(FatalError);
        }

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const faceList& faces = mesh_.faces();

        auto tmeshPhi(const_cast<fvMesh&>(mesh_).setPhi());
        if (tmeshPhi)
        {
            const scalar rdt = 1.0/mesh_.time().deltaTValue();

            // Set the mesh flux
            auto& meshPhi = tmeshPhi.ref();
            auto& meshPhii = meshPhi.primitiveFieldRef();
            auto& meshPhiBf = meshPhi.boundaryFieldRef();

            //meshPhi == dimensionedScalar(dimVolume/dimTime, Zero);
            meshPhii = Zero;
            meshPhiBf == Zero;

            forAll(changedFaceIDs_, i)
            {
                const face& f = faces[changedFaceIDs_[i]];

                if (changedPatchIDs_[i] == -1)
                {
                    const label facei = changedFaceIDs_[i];
                    meshPhii[facei] = f.sweptVol(oldPoints, currPoints)*rdt;
                }
                else
                {
                    const label patchi = changedPatchIDs_[i];
                    const polyPatch& pp = pbm[patchi];

                    if (isA<emptyPolyPatch>(pp))
                    {
                        continue;
                    }

                    const label patchFacei = changedFaceIDs_[i] - pp.start();

                    meshPhiBf[patchi][patchFacei] =
                        f.sweptVol(oldPoints, currPoints)*rdt;
                }
            }
        }

        if (partialUpdate_ && haveGeometry)
        {
            // Keep base geometry and update as needed
            DebugInFunction << "Performing partial geometry update" << endl;

            // Initialise geometry using the old/existing values
            vectorField faceCentres(mesh_.faceCentres());
            vectorField faceAreas(mesh_.faceAreas());

            // Make face centres and areas consistent with new points
            primitiveMeshTools::updateFaceCentresAndAreas
            (
                mesh_,
                changedFaceIDs_,
                mesh_.points(),
                faceCentres,
                faceAreas
            );

            vectorField cellCentres(mesh_.cellCentres());
            scalarField cellVolumes(mesh_.cellVolumes());

            primitiveMeshTools::updateCellCentresAndVols
            (
                mesh_,
                faceCentres,
                faceAreas,
                changedCellIDs_,
                mesh_.cells(),
                cellCentres,
                cellVolumes
            );

            const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
            (
                std::move(faceCentres),
                std::move(faceAreas),
                std::move(cellCentres),
                std::move(cellVolumes)
            );

            if (debug)
            {
                for (const auto& p : mesh_.boundaryMesh())
                {
                    Pout<< "SBM --- " << p.name()
                        << " sum(Sf)=" << sum(p.faceAreas())
                        << " sum(mag(Sf))=" << sum(mag(p.faceAreas()))
                        << endl;
                }
            }
        }
        else
        {
            DebugInFunction
                << "Performing complete geometry clear and update" << endl;

            // Clear out old geometry
            // Note: this recreates the old primitiveMesh::movePoints behaviour
            const_cast<fvMesh&>(mesh_).primitiveMesh::clearGeom();

            // Use lower level to calculate the geometry
            const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
        }
    }
    else
    {
        DebugInFunction << "Performing complete geometry update" << endl;

        // Use lower level to calculate the geometry
        const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
    }
}


void Foam::solidBodyFvGeometryScheme::updateMesh(const mapPolyMesh& mpm)
{
    cacheInitialised_ = false;
}


bool Foam::solidBodyFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    if
    (
        (faceCentres.size() != mesh_.nFaces())
     || (faceAreas.size() != mesh_.nFaces())
     || (cellCentres.size() != mesh_.nCells())
     || (cellVolumes.size() != mesh_.nCells())
     || !oldPoints
    )
    {
        // Do all
        return basicFvGeometryScheme::updateGeom
        (
            points,
            oldPoints,
            faceCentres,
            faceAreas,
            cellCentres,
            cellVolumes
        );
    }
    else
    {
        // Since oldPoints provided assume that face & cell geometry is
        // up to date with it

        bitSet isChangedPoint;
        bitSet isChangedFace;
        bitSet isChangedCell;
        const bool changed = markChanges
        (
            oldPoints(),
            points,
            isChangedPoint,
            isChangedFace,
            isChangedCell
        );

        if (!changed)
        {
            return false;
        }

        // Make face centres and areas consistent with new points
        primitiveMeshTools::updateFaceCentresAndAreas
        (
            mesh_,
            isChangedFace.toc(),
            points,
            faceCentres,
            faceAreas
        );

        primitiveMeshTools::updateCellCentresAndVols
        (
            mesh_,
            faceCentres,
            faceAreas,
            isChangedCell.toc(),
            mesh_.cells(),
            cellCentres,
            cellVolumes
        );

        return true;
    }
}


// ************************************************************************* //
