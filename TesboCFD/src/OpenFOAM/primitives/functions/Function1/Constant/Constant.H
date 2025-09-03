/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::Function1Types::Constant

Description
    Templated function that returns a constant value.

    Usage - for entry \<entryName\> returning the value \<value\>,
    can be specified is different formats.

    Inline specification:
    \verbatim
        <entryName>  constant <value>
    \endverbatim

    Dictionary format:
    \verbatim
        <entryName>
        {
            type    constant;
            value   <value>;
        }
    \endverbatim

SourceFiles
    Constant.C
    ConstantI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Constant_H
#define Foam_Function1Types_Constant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public Function1<Type>
{
    // Private Data

        //- Constant value
        Type value_;


public:

    // Runtime type information
    TypeName("constant");


    // Generated Methods

        //- No copy assignment
        void operator=(const Constant<Type>&) = delete;


    // Constructors

        //- Construct from components
        Constant
        (
            const word& entryName,
            const Type& value,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from entry name, dictionary and optional registry
        Constant
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from entry name and Istream
        //  Reads the constant value without the Function1 type
        //  for backward compatibility
        Constant(const word& entryName, Istream& is);

        //- Copy constructor
        explicit Constant(const Constant<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Constant() = default;


    // Member Functions

        //- Change the constant value
        inline void reset(const Type& val);

        //- Value is independent of x
        virtual inline bool constant() const { return true; }

        //- Return constant value
        virtual inline Type value(const scalar) const;

        //- Integrate between two values
        virtual inline Type integrate(const scalar x1, const scalar x2) const;

        //- Return value as a function of (scalar) independent variable
        virtual tmp<Field<Type>> value(const scalarField& x) const;

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ConstantI.H"

#ifdef NoRepository
    #include "Constant.C"
    #include "Function1New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
