/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::solidBodyMotionFunctions::oscillatingLinearMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating displacement.

    The governing equation for the displacement is given by:
    \f[
        y =  A*\sin(B*(t + C)) + D
    \f]

    where:
    \vartable
        y    | Displacement [m] (vector)
        A    | Amplitude [m] (vector)
        B    | Radial velocity [rad/s] (scalar)
        C    | Phase-shift to left [s] (scalar)
        D    | Vertical shift [m] (vector)
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    oscillation
    {
        // Mandatory entries
        solidBodyMotionFunction oscillatingLinearMotion;

        oscillatingLinearMotionCoeffs
        {
            // Mandatory entries
            amplitude       <Function1<vector>>;
            omega           <Function1<scalar>>;

            // Optional entries
            phaseShift      <Function1<scalar>>;
            verticalShift   <Function1<vector>>;
        }
    }

    where the entries mean:
    \table
      Property  | Description          | Type                | Reqd | Deflt
      solidBodyMotionFunction          | Type name: oscillatingLinearMotion <!--
                -->                    | word                | yes  | -
      amplitude | Amplitude            | Function1\<vector\> | yes  | -
      omega     | Radial velocity      | Function1\<scalar\> | yes  | -
      phaseShift | Phase-shift to left | Function1\<scalar\> | no   | -
      verticalShift | Vertical shift   | Function1\<vector\> | no   | -
    \endtable

SourceFiles
    oscillatingLinearMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidBodyMotionFunctions_oscillatingLinearMotion_H
#define Foam_solidBodyMotionFunctions_oscillatingLinearMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                    Class oscillatingLinearMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingLinearMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Radial velocity
        autoPtr<Function1<scalar>> omegaPtr_;

        //- Phase shift
        autoPtr<Function1<scalar>> phaseShiftPtr_;

        //- Amplitude
        autoPtr<Function1<vector>> amplitudePtr_;

        //- Vertical shift
        autoPtr<Function1<vector>> verticalShiftPtr_;


    // Private Member Functions

        //- No copy construct
        oscillatingLinearMotion(const oscillatingLinearMotion&) = delete;

        //- No copy assignment
        void operator=(const oscillatingLinearMotion&) = delete;


public:

    //- Runtime type information
    TypeName("oscillatingLinearMotion");


    // Constructors

        //- Construct from components
        oscillatingLinearMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new oscillatingLinearMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~oscillatingLinearMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
