/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::solidBodyMotionFunctions::drivenLinearMotion

Description
    Variable velocity displacement. The velocity is read from a
    uniformVectorField from the time registry with the name CofGvelocity.

    Optional plane of motion can be added with the normal vector


SourceFiles
    drivenLinearMotion.C

\*---------------------------------------------------------------------------*/

#ifndef drivenLinearMotion_H
#define drivenLinearMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class drivenLinearMotion Declaration
\*---------------------------------------------------------------------------*/

class drivenLinearMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Name of the meshObject to dum CofG displacement
        word cOfGdisplacement_;

        //- Uniform vector to follow
        uniformDimensionedVectorField CofGdisp_;


    // Private Member Functions

        //- No copy construct
        drivenLinearMotion(const drivenLinearMotion&) = delete;

        //- No copy assignment
        void operator=(const drivenLinearMotion&) = delete;


public:

    //- Runtime type information
    TypeName("drivenLinearMotion");


    // Constructors

        //- Construct from components
        drivenLinearMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new drivenLinearMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~drivenLinearMotion() = default;


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
