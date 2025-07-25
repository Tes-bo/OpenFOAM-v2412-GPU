/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::sixDoFRigidBodyMotionConstraints::axis

Description
    This constraint imposes an orientation limitation where bodies
    are restricted to rotate only around a fixed axis, optionally with
    the inclusion of maximum and minimum rotation angle constraints.

Usage
    Minimal example by \c constant/dynamicMeshDict.sixDoFRigidBodyMotionCoeffs:
    \verbatim
    constraints
    {
        constrainRotation1
        {
            // Mandatory entries
            sixDoFRigidBodyMotionConstraint     axis;
            axis                                <vector>;

            // Optional entries
            maxClockwiseTheta                   <Function1<scalar>>;
            maxCounterclockwiseTheta            <Function1<scalar>>;
            thetaUnits                          <word>;
            referenceOrientation                <tensor>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      sixDoFRigidBodyMotionConstraint | Type name: axis | word | yes  | -
      axis         | Reference rotation axis fixed in global space   <!--
                   -->                                | vector | yes  | -
      maxClockwiseTheta | Maximum clockwise rotation angle about the <!--
                   --> rotation axis            | Function1\<scalar\> | no | -
      maxCounterlockwiseTheta | Maximum counterclockwise rotation    <!--
             --> angle about the rotation axis  | Function1\<scalar\> | no | -
      thetaUnits   | Units of theta angles             | word   | no  | radians
      referenceOrientation | Reference orientation where there is no moment <!--
                   -->                                 | tensor | no  | I
    \endtable

    The inherited entries are elaborated in:
      - \link Function1.H \endlink

Notes
  - The units for \c thetaUnits can be specified
    as either \c degrees or \c radians.
  - The \c maxClockwiseTheta and \c maxCounterlockwiseTheta
    are always non-negative.
  - Negative and positive \c theta correspond to clockwise
    and counterclockwise rotation sectors with respect to
    the reference orientation, respectively.

SourceFiles
    sixDoFRigidBodyMotionAxisConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sixDoFRigidBodyMotionAxisConstraint_H
#define Foam_sixDoFRigidBodyMotionAxisConstraint_H

#include "sixDoFRigidBodyMotionConstraint.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionConstraints
{

/*---------------------------------------------------------------------------*\
                            Class axis Declaration
\*---------------------------------------------------------------------------*/

class axis
:
    public sixDoFRigidBodyMotionConstraint
{
    // Private Data

        //- Reference orientation where there is no moment
        tensor refQ_;

        //- Reference rotation axis fixed in global space
        vector axis_;

        //- Maximum clockwise rotation angle around the rotation axis
        autoPtr<Function1<scalar>> maxCWThetaPtr_;

        //- Maximum counterclockwise rotation angle around the rotation axis
        autoPtr<Function1<scalar>> maxCCWThetaPtr_;

        //- Flag to enable to use maxTheta in degrees rather than in radians
        bool degrees_;


    // Private Member Functions

        //- Return the index of the rotation sector that the
        //- body resides with respect to the reference orientation
        //  -1: clockwise, 0: reference orientation, 1: counterclockwise
        label rotationSector(const vector& oldDir, const vector& newDir) const;

        //- Return the direction of the given force or momentum
        //  false: clockwise, true: counterclockwise
        bool calcDir(const vector& fm, const bool rotationSector) const;


public:

    //- Runtime type information
    TypeName("axis");


    // Constructors

        //- Construct from components
        axis
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const sixDoFRigidBodyMotion& motion
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionConstraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionConstraint>
            (
                new axis(*this)
            );
        }

        //- No copy assignment
        void operator=(const sixDoFRigidBodyMotionConstraint&) = delete;


    //- Destructor
    virtual ~axis() = default;


    // Member Functions

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
