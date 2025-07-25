/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Foam::Function1Types::Scale

Description
    Function1 which scales a given 'value' function by a scalar 'scale'
    function.

    This is particularly useful to ramp a time-varying value by one of the
    monotonic ramp functions.

    Usage for a vector:
    \verbatim
        <entryName>
        {
            type      scale;

            scale
            {
                type        linearRamp;

                start       0;
                duration    10;
            }

            value
            {
                type        sine;

                frequency   10;
                amplitude   1;
                scale       (1 0.1 0);
                level       (10 1 0);
            }
        }
    \endverbatim

    Where:
    \table
        Property | Description                                  | Required
        scale    | Scaling function of type Function1<scalar>   | yes
        value    | Function of type Function1<Type>             | yes
    \endtable

SourceFiles
    Scale.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Scale_H
#define Foam_Function1Types_Scale_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Scale Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Scale
:
    public Function1<Type>
{
    // Private Data

        //- Scalar scaling function
        autoPtr<Function1<scalar>> scale_;

        //- Value function
        autoPtr<Function1<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    //- Runtime type information
    TypeName("scale");


    // Generated Methods

        //- No copy assignment
        void operator=(const Scale<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Scale
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit Scale(const Scale<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Scale() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ScaleI.H"

#ifdef NoRepository
    #include "Scale.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
