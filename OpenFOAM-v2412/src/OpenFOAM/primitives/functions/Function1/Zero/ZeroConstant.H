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
    Foam::Function1Types::ZeroConstant

Description
    Templated function that returns the corresponding 0 (zero).

    Usage:
    Inline specification:
    \verbatim
        <entryName>  zero;
    \endverbatim

    In dictionary format:
    \verbatim
        <entryName>
        {
            type    zero;
        }
    \endverbatim

SourceFiles
    ZeroConstant.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_ZeroConstant_H
#define Foam_Function1Types_ZeroConstant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class ZeroConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ZeroConstant
:
    public Function1<Type>
{
public:

    //- Runtime type information
    TypeName("zero");


    // Generated Methods

        //- Default copy construct
        ZeroConstant(const ZeroConstant<Type>&) = default;

        //- No copy assignment
        void operator=(const ZeroConstant<Type>&) = delete;


    // Constructors

        //- Construct from entry name and optional registry
        explicit ZeroConstant
        (
            const word& entryName,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from entry name, dictionary and optional registry
        ZeroConstant
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
    virtual ~ZeroConstant() = default;


    // Member Functions

        //- Value is independent of x
        virtual inline bool constant() const { return true; }

        //- Return constant value
        virtual inline Type value(const scalar) const
        {
            return pTraits<Type>::zero;
        }

        //- Integrate between two values
        virtual inline Type integrate(const scalar x1, const scalar x2) const
        {
            return pTraits<Type>::zero;
        }

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ZeroConstant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
