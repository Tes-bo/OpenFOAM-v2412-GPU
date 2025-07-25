/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::Function1Types::Square

Description
    A templated square-wave function with support for offset, etc.

    The wave period can be specified directly

    \f[
        a square((t - t0) / p)) s + l
    \f]

    Or it can be specified by the frequency

    \f[
        a square(f (t - t0)) s + l
    \f]

    where

    \f$ square(t) \f$ is the square-wave function in range \f$ [-1, 1] \f$
    with a mark/space ratio of \f$ r \f$

    \vartable
        symbol  | Description       | Units
        a       | Amplitude         | -
        f       | Frequency         | [1/s]
        p       | Period            | [s]
        s       | Type scale factor | -
        l       | Type offset level | -
        t       | Time              | [s]
        t0      | Start time offset | [s]
        r       | mark/space ratio  | -
    \endvartable

    The dictionary specification would typically resemble this:
    \verbatim
    entry1
    {
        type        square;
        frequency   10;
        amplitude   0.1;

        // A scalar Function1
        scale       2e-6;
        level       2e-6;
    }
    entry2
    {
        type        square;
        frequency   10;

        // A vector Function1
        scale       (1 0.1 0);
        level       (10 1 0);
    }
    \endverbatim

    where the entries mean:
    \table
        Property  | Description             | Type | Reqd | Default
        type      | Function type: square   | word | yes |
        amplitude | Amplitude               | Function1<scalar> | no | 1
        frequency | Frequency [1/s]         | Function1<scalar> | or period |
        period    | Period [s]              | Function1<scalar> | or frequency |
        scale     | Scale factor (Type)     | Function1<Type> | yes |
        level     | Offset level (Type)     | Function1<Type> | yes |
        t0        | Start time offset       | scalar       | no | 0
        mark      | Positive amount         | scalar       | no | 1
        space     | Negative amount         | scalar       | no | 1
    \endtable

Note
    For slow oscillations it can be more intuitive to specify the period.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Square_H
#define Foam_Function1Types_Square_H

#include "Sine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Square Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Square
:
    public Function1Types::Sine<Type>
{
    // Private Data

        //- Mark/space ratio (square function only)
        scalar mark_;

        //- Mark/space ratio (square function only)
        scalar space_;

public:

    // Runtime type information
    TypeName("square");


    // Generated Methods

        //- No copy assignment
        void operator=(const Square<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Square
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit Square(const Square<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Square() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const
        {
            return Sine<Type>::squareValue(t, mark_ / (mark_ + space_));
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Square.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
