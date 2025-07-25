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
    Foam::Function1Types::OneConstant

Description
    Templated function that returns the corresponding 1 (one).

    Usage:
    \verbatim
        <entryName> one;
    \endverbatim

SourceFiles
    OneConstant.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_OneConstant_H
#define Foam_Function1Types_OneConstant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class OneConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class OneConstant
:
    public Function1<Type>
{
public:

    //- Runtime type information
    TypeName("one");


    // Generated Methods

        //- Default copy construct
        OneConstant(const OneConstant<Type>&) = default;

        //- No copy assignment
        void operator=(const OneConstant<Type>&) = delete;


    // Constructors

        //- Construct from entry name and optional registry
        explicit OneConstant
        (
            const word& entryName,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from entry name, dictionary and optional registry
        OneConstant
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~OneConstant() = default;


    // Member Functions

        //- Value is independent of x
        virtual inline bool constant() const { return true; }

        //- Return constant value
        virtual inline Type value(const scalar) const
        {
            return pTraits<Type>::one;
        }

        //- Integrate between two values
        virtual inline Type integrate(const scalar x1, const scalar x2) const
        {
            return (x2 - x1)*pTraits<Type>::one;
        }

        //- Return value as a function of (scalar) independent variable
        virtual tmp<Field<Type>> value(const scalarField& x) const;

        //- Integrate between two (scalar) values
        virtual tmp<Field<Type>> integrate
        (
            const scalarField& x1,
            const scalarField& x2
        ) const;

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OneConstant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
