/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotionAxisConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(axis, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        axis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sixDoFRigidBodyMotionConstraints::axis::rotationSector
(
    const vector& oldDir,
    const vector& newDir
) const
{
    const scalar thetaDir = (oldDir ^ newDir) & axis_;

    if (equal(thetaDir, 0))
    {
        return 0;
    }

    return label(sign(thetaDir));
}


bool Foam::sixDoFRigidBodyMotionConstraints::axis::calcDir
(
    const vector& fm,
    const bool rotationSector
) const
{
    const scalar fmDir = axis_ & fm;

    if (equal(fmDir, 0))
    {
        return rotationSector;
    }

    return (label(sign(fmDir)) == 1) ? true : false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::axis::axis
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotion& motion
)
:
    sixDoFRigidBodyMotionConstraint(name, sDoFRBMCDict, motion),
    refQ_(),
    axis_(),
    maxCWThetaPtr_(),
    maxCCWThetaPtr_(),
    degrees_(false)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionConstraints::axis::constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::sixDoFRigidBodyMotionConstraints::axis::constrainRotation
(
    pointConstraint& pc
) const
{
    if (!(maxCWThetaPtr_ && maxCCWThetaPtr_))
    {
        pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
        return;
    }


    // Calculate principal directions of the body
    const vector refDir
    (
        rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 1, 0)
    );
    const vector oldDir
    (
        (refQ_ & refDir).removeCollinear(axis_).normalise()
    );
    const vector newDir
    (
        (motion().orientation() & refDir).removeCollinear(axis_).normalise()
    );


    // Find the index of the rotation sector that the body resides
    const label rotationSectorIndex = rotationSector(oldDir, newDir);

    if (!rotationSectorIndex)
    {
        // The body resides at the reference orientation
        pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
        return;
    }

    const bool rotationSector = (rotationSectorIndex == 1) ? true : false;


    // Calculate the directions of total momentum and force acting on the body
    const bool angularMomentumDir =
        calcDir
        (
            motion().state().pi(),
            rotationSector
        );
    const bool torqueDir = calcDir(motion().state().tau(), rotationSector);


    // Calculate the rotation angle of the body wrt the reference orientation
    const scalar theta = mag(acos(min(oldDir & newDir, scalar(1))));


    // Calculate maximum clockwise and counterclockwise rotation angles
    const scalar t = motion().time().timeOutputValue();
    const scalar maxCWTheta =
        degrees_
      ? mag(degToRad(maxCWThetaPtr_->value(t)))
      : mag(maxCWThetaPtr_->value(t));

    const scalar maxCCWTheta =
        degrees_
      ? mag(degToRad(maxCCWThetaPtr_->value(t)))
      : mag(maxCCWThetaPtr_->value(t));


    // Apply the constraints according to various conditions
    if
    (
        (rotationSector && (theta < maxCCWTheta))
     || (!rotationSector && (theta < maxCWTheta))
    )
    {
        pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
    }
    else
    {
        if (rotationSector == angularMomentumDir)
        {
            const scalar magPi = mag(motion().state().pi());

            if (equal(magPi, scalar(0)) && rotationSector != torqueDir)
            {
                pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
            }
            else
            {
                // Constrain all rotational motions
                pc.combine(pointConstraint(Tuple2<label, vector>(3, Zero)));
            }
        }
        else
        {
            // If there is a difference between the directions of
            // body rotation and of torque, release the body
            pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
        }
    }


    if (motion().report())
    {
        Info
            << " old direction = " << oldDir << nl
            << " new direction = " << newDir << nl
            << " rotationSector = " << rotationSector << nl
            << " theta = " << sign((oldDir ^ newDir) & axis_)*theta << nl
            << " torque = " << motion().state().tau() << nl
            << " torque dir = " << torqueDir << nl
            << " angular momentum = " << motion().state().pi() << nl
            << " angular momentum dir = " << angularMomentumDir << nl
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionConstraints::axis::read
(
    const dictionary& sDoFRBMCDict
)
{
    if (!sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict))
    {
        return false;
    }

    sDoFRBMCCoeffs_.readEntry("axis", axis_);

    axis_.normalise();

    if (mag(axis_) < VSMALL)
    {
        FatalIOErrorInFunction(sDoFRBMCDict)
            << "The entry 'axis' cannot have zero length: " << axis_
            << exit(FatalIOError);
    }


    if (sDoFRBMCCoeffs_.found("thetaUnits"))
    {
        const word thetaUnits(sDoFRBMCCoeffs_.getWord("thetaUnits"));

        if (thetaUnits == "degrees")
        {
            degrees_ = true;
        }
        else if (thetaUnits == "radians")
        {
            degrees_ = false;
        }
        else
        {
            FatalIOErrorInFunction(sDoFRBMCCoeffs_)
                << "The units of thetaUnits can be either degrees or radians"
                << exit(FatalIOError);
        }


        maxCWThetaPtr_.reset
        (
            Function1<scalar>::New
            (
                "maxClockwiseTheta",
                sDoFRBMCCoeffs_,
                &motion().time()
            )
        );

        maxCCWThetaPtr_.reset
        (
            Function1<scalar>::New
            (
                "maxCounterclockwiseTheta",
                sDoFRBMCCoeffs_,
                &motion().time()
            )
        );


        refQ_ =
            sDoFRBMCCoeffs_.getOrDefault<tensor>("referenceOrientation", I);

        if (mag(mag(refQ_) - sqrt(3.0)) > ROOTSMALL)
        {
            FatalIOErrorInFunction(sDoFRBMCCoeffs_)
                << "The entry 'referenceOrientation' " << refQ_
                << " is not a rotation tensor. "
                << "mag(referenceOrientation) - sqrt(3) = "
                << mag(refQ_) - sqrt(3.0) << nl
                << exit(FatalIOError);
        }
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::axis::write
(
    Ostream& os
) const
{
    os.writeEntry("axis", axis_);


    if (maxCWThetaPtr_ && maxCCWThetaPtr_)
    {
        if (degrees_)
        {
            os.writeEntry("thetaUnits", "degrees");
        }
        else
        {
            os.writeEntry("thetaUnits", "radians");
        }

        if (maxCWThetaPtr_)
        {
            maxCWThetaPtr_->writeData(os);
        }

        if (maxCCWThetaPtr_)
        {
            maxCCWThetaPtr_->writeData(os);
        }

        os.writeEntry("referenceOrientation", refQ_);
    }
}


// ************************************************************************* //
