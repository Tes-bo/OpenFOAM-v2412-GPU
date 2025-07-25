/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "lumpedPointDisplacementPointPatchVectorField.H"
#include "lumpedPointMovement.H"
#include "lumpedPointIOMovement.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "Time.H"
#include "polyMesh.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        lumpedPointDisplacementPointPatchVectorField
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::lumpedPointDisplacementPointPatchVectorField::setPatchControls
(
    const pointVectorField& pvf,
    const pointField& points0
)
{
    label count = 0;

    const pointVectorField::Boundary& bf = pvf.boundaryField();
    const polyBoundaryMesh& patches = pvf.mesh().mesh().boundaryMesh();

    forAll(bf, patchi)
    {
        // Patch of this type
        const auto* p = isA<patchType>(bf[patchi]);

        if (p)
        {
            // Patch controls (mapping) for calculating forces/moments
            const_cast<lumpedPointMovement&>(p->movement())
                .setPatchControl
                (
                    patches[patchi],
                    p->controllers(),
                    points0
                );

            ++count;
        }
    }

    return count;
}


Foam::label Foam::lumpedPointDisplacementPointPatchVectorField::setInterpolators
(
    const pointVectorField& pvf,
    const pointField& points0
)
{
    label count = 0;

    const pointVectorField::Boundary& bf = pvf.boundaryField();

    forAll(bf, patchi)
    {
        // Patch of this type
        const auto* p = isA<patchType>(bf[patchi]);

        if (p)
        {
            // Patch controls (mapping) for calculating forces/moments
            const_cast<lumpedPointMovement&>(p->movement())
                .setInterpolator
                (
                    p->patch(),
                    points0
                );

            ++count;
        }
    }

    return count;
}


Foam::labelList
Foam::lumpedPointDisplacementPointPatchVectorField::patchIds
(
    const pointVectorField& pvf
)
{
    const pointVectorField::Boundary& bf = pvf.boundaryField();

    DynamicList<label> patchLst(bf.size());
    forAll(bf, patchi)
    {
        // Patch of this type
        if (isA<patchType>(bf[patchi]))
        {
            patchLst.append(patchi);
            // or patchLst.append(bf[patchi].patch().index());
        }
    }

    return patchLst.shrink();
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::pointField&
Foam::lumpedPointDisplacementPointPatchVectorField::points0() const
{
    const objectRegistry& obr = this->patch().boundaryMesh().mesh().db();

    // Obtain starting locations from the motionSolver (when possible)
    const auto* solver =
        obr.cfindObject<displacementMotionSolver>("dynamicMeshDict");

    if (solver)
    {
        if (points0Ptr_)
        {
            points0Ptr_.reset(nullptr);
        }
        return solver->points0();
    }
    else if (!points0Ptr_)
    {
        points0Ptr_.reset
        (
            new pointIOField
            (
                points0MotionSolver::points0IO
                (
                    this->patch().boundaryMesh().mesh().mesh()
                )
            )
        );
    }

    return *points0Ptr_;
}


const Foam::lumpedPointMovement&
Foam::lumpedPointDisplacementPointPatchVectorField::movement() const
{
    const objectRegistry& obr = this->patch().boundaryMesh().mesh().db();

    lumpedPointIOMovement* ptr =
        lumpedPointIOMovement::getMovementObject(obr);

    if (!ptr)
    {
        // Create and register with this patch as the owner
        ptr = lumpedPointIOMovement::New(obr, this->patch().index()).ptr();

        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    controllers_(),
    dataWritten_(0, 0),
    points0Ptr_(nullptr)
{}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    controllers_(),
    dataWritten_(0, 0),
    points0Ptr_(nullptr)
{
    dict.readIfPresent("controllers", controllers_);

    dict.readIfPresent("dataWritten", dataWritten_);

    if (controllers_.empty())
    {
        WarningInFunction
            << "No controllers specified, using all lumped points for patch: "
            << this->patch().name() << nl << nl;
    }

    // controllers_ : check? remove duplicates?
}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const lumpedPointDisplacementPointPatchVectorField& rhs,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(rhs, p, iF, mapper),
    controllers_(rhs.controllers_),
    dataWritten_(rhs.dataWritten_),
    points0Ptr_(nullptr)
{}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const lumpedPointDisplacementPointPatchVectorField& rhs,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(rhs, iF),
    controllers_(rhs.controllers_),
    dataWritten_(rhs.dataWritten_),
    points0Ptr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lumpedPointDisplacementPointPatchVectorField::
~lumpedPointDisplacementPointPatchVectorField()
{
    // de-register movement if in use and managed by this patch
    lumpedPointIOMovement* ptr =
        lumpedPointIOMovement::getMovementObject
        (
            this->patch().boundaryMesh().mesh().db()
        );

    if (ptr && ptr->ownerId() == this->patch().index())
    {
        movement().coupler().shutdown();

        ptr->checkOut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label timeIndex = this->db().time().timeIndex();

    enum Time::stopAtControls action = Time::stopAtControls::saUnknown;

    if (movement().ownerId() == this->patch().index())
    {
        // The ownerId is always the lowest patch number,
        // thus will always be triggered first

        if (lumpedPointMovement::debug)
        {
            Pout<<"masterPatch: " << this->patch().index() << endl;
        }

        const polyMesh& mesh = this->patch().boundaryMesh().mesh().mesh();

        // Mapping for calculating forces
        // as well as face point interpolation
        if (!movement().hasMapping())
        {
            // Add mapping for calculating forces/moments
            setPatchControls
            (
                static_cast<const pointVectorField&>
                (
                    this->internalField()
                ),
                this->points0()
            );
        }

        int triggered = 0;

        if
        (
            movement().coupler().slaveFirst()
         && !movement().coupler().initialized()
        )
        {
            // Master does nothing yet, but slave will
            triggered = -1;
        }
        else if (movement().couplingPending(timeIndex))
        {
            // Trigger is pending, or coupling not yet initialized
            triggered = 1;
        }

        if (triggered > 0)
        {
            // Synchronized for all processes
            List<vector> forces, moments;
            movement().forcesAndMoments(mesh, forces, moments);

            if (lumpedPointMovement::debug)
            {
                Pout<<"gatherForces: " << forces << " called from patch "
                    << this->patch().index() << endl;

                Info<< "output forces to file: called from patch "
                    << this->patch().index() << nl
                    << "# " << forces.size() << " force entries" << nl
                    << "# fx fy fz" << nl
                    << "output forces to file: "
                    << forces << " called from patch "
                    << this->patch().index() << endl;
            }

            // Update times when data (forces) were written
            // With first=time, second=prevTime

            dataWritten_.second() = dataWritten_.first();
            dataWritten_.first() = this->db().time().timeOutputValue();

            if (Pstream::master())
            {
                movement().writeData(forces, moments, &dataWritten_);

                // Signal external source to execute
                movement().coupler().useSlave();
            }
        }

        if (triggered)
        {
            // Wait for slave to provide data (includes MPI barrier)
            // and catch any abort information sent from slave
            action = movement().coupler().waitForSlave();

            const_cast<lumpedPointMovement&>(movement()).readState();

            movement().couplingCompleted(timeIndex);
        }
    }

    if (!movement().hasInterpolator(this->patch()))
    {
        const_cast<lumpedPointMovement&>(movement()).setInterpolator
        (
            this->patch(),
            this->points0()
        );
    }

    tmp<pointField> tdisp =
        movement().pointsDisplacement
        (
            this->patch(),
            this->points0()
        );

    this->operator==(tdisp);


    fixedValuePointPatchField<vector>::updateCoeffs();

    // Process any abort information sent from slave
    if
    (
        action != this->db().time().stopAt()
     && action != Time::stopAtControls::saUnknown
    )
    {
        this->db().time().stopAt(action);
    }
}


void Foam::lumpedPointDisplacementPointPatchVectorField::write(Ostream& os)
const
{
    pointPatchField<vector>::write(os);

    if (controllers_.size())
    {
        os.writeEntry("controllers", controllers_);
    }

    // Times when data were written is only meaningful on the owner patch
    if (movement().ownerId() == this->patch().index())
    {
        os.writeEntry("dataWritten", dataWritten_);
    }

    this->writeValueEntry(os);
}


// ************************************************************************* //
