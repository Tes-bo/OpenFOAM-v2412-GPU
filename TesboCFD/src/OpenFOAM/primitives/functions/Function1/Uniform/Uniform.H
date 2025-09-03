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
    Foam::Function1Types::Uniform

Description
    Templated function that returns a constant value.

    Provides backward-compatibility for cases where a field is spatially
    "uniform" and may be treated as a constant value.

    Usage - for entry \<entryName\> returning the value <value>:
    \verbatim
        <entryName>    uniform  <value>
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Uniform_H
#define Foam_Function1Types_Uniform_H

#include "Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Uniform Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Uniform
:
    public Constant<Type>
{
public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("uniform");


    // Generated Methods

        //- Default copy construct
        Uniform(const Uniform<Type>&) = default;

        //- No copy assignment
        void operator=(const Uniform<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Uniform
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        )
        :
            Constant<Type>(entryName, dict, obrPtr)
        {}

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
