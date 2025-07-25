/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "dynamicMotionSolverFvMeshAMI.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "polyTopoChange.H"
#include "MeshObject.H"
#include "lduMesh.H"
#include "surfaceInterpolate.H"

#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverFvMeshAMI, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMeshAMI,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMeshAMI,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverFvMeshAMI::dynamicMotionSolverFvMeshAMI
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMotionSolverFvMeshAMI::init(const bool doInit)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    motionPtr_ = motionSolver::New(*this);

    // allow restarts during initialization to match patch field values if
    // required 
    const auto& pbm = boundaryMesh();
    bool changeRequired = false;
    for (label patchi = 0; patchi < pbm.nNonProcessor(); ++patchi)
    {
        const auto* cycAmiPtr = isA<cyclicAMIPolyPatch>(pbm[patchi]);

        if (cycAmiPtr)
        {
            changeRequired = cycAmiPtr->createAMIFaces() || changeRequired;
        }
        else
        {
            const auto* cycAcmiPtr = isA<cyclicACMIPolyPatch>(pbm[patchi]);
            if (cycAcmiPtr)
            {
                changeRequired =
                    cycAcmiPtr->cyclicAMIPolyPatch::createAMIFaces()
                 || changeRequired;
            }
        }
    }
    if (returnReduceOr(changeRequired))
    {
        update();
    }

    return true;
}


Foam::dynamicMotionSolverFvMeshAMI::dynamicMotionSolverFvMeshAMI
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    motionPtr_(motionSolver::New(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::dynamicMotionSolverFvMeshAMI::motion() const
{
    return *motionPtr_;
}


bool Foam::dynamicMotionSolverFvMeshAMI::update()
{
    // Mesh not moved/changed yet
    moving(false);
    topoChanging(false);

    if (debug)
    {
        for (const fvPatch& fvp : boundary())
        {
            if (!isA<processorFvPatch>(fvp))
            {
                Info<< "1 --- patch:" << fvp.patch().name()
                    << " area:" << gSum(fvp.magSf()) << endl;
            }
        }
    }

    pointField newPoints(motionPtr_->curPoints());

    polyBoundaryMesh& pbm = const_cast<polyBoundaryMesh&>(boundaryMesh());

    // Scan all patches and see if we want to apply a mesh topology  update
    bool changeRequired = false;
    for (label patchi = 0; patchi < pbm.nNonProcessor(); ++patchi)
    {
        const polyPatch& pp = pbm[patchi];

        DebugInfo
            << "pre-topology change: patch " << pp.name()
            << " size:" << returnReduce(pp.size(), sumOp<label>())
            << " mag(faceAreas):" << gSum(mag(pp.faceAreas())) << endl;

        //changeRequired = pp.changeTopology(newPoints) || changeRequired;
        changeRequired = pp.changeTopology() || changeRequired;
    }

    if (returnReduceOr(changeRequired))
    {
        polyTopoChange polyTopo(*this);

        // Set new point positions in polyTopo object
        polyTopo.movePoints(newPoints);

        // Accumulate the patch-based mesh changes on the current mesh
        // Note:
        // - updates the AMIs using the new points
        // - creates a topo change object that removes old added faces and
        //   adds the new faces
        for (polyPatch& pp : pbm)
        {
            pp.setTopology(polyTopo);
        }

        // Update geometry
        // Note
        // - changeMesh leads to polyMesh::resetPrimitives which will also
        //   trigger polyBoundaryMesh::updateMesh (init and update) and
        //   ::calcGeometry (with topoChanging = false)
        // - BUT: mesh still corresponds to original (non-extended mesh) so
        //   we want to bypass these calls...
        // - after changes topoChanging = true
        autoPtr<mapPolyMesh> map =
            polyTopo.changeMesh
            (
                *this,
                true       // We will be calling movePoints after this update
            );

        // Apply topology change - update fv geometry and map fields
        // - polyMesh::updateMesh
        //   - fires initUpdateMesh and updateMesh in AMI BCs - called before
        //     mapFields
        // - AMI addressing must be up-to-date - used by, e.g. FaceCellWave
        // - will trigger (again) polyBoundaryMesh::updateMesh (init and update)
        updateMesh(map());

        // Move points and update derived properties
        // Note:
        // - resets face areas based on raw point locations!
        // - polyBoundaryMesh::updateMesh (init and update)
        // Note:
        // - processorPolyPatches will trigger calculation of faceCentres
        //   (and therefore cell volumes), so need to update faceAreas in
        //   initMovePoints since proc patches will be evaluated later than
        //   AMI patches
        if (map().hasMotionPoints())
        {
            movePoints(map().preMotionPoints());
        }
    }
    else
    {
        fvMesh::movePoints(newPoints);
    }

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();

        surfaceVectorField* UfPtr = getObjectPtr<surfaceVectorField>("Uf");
        if (UfPtr)
        {
            *UfPtr = fvc::interpolate(*Uptr);
        }
    }

    if (debug)
    {
        for (const fvPatch& fvp : boundary())
        {
            if (!isA<processorFvPatch>(fvp))
            {
                Info<< "2 --- patch:" << fvp.patch().name()
                    << " area:" << gSum(fvp.magSf()) << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
