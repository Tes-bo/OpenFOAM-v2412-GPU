/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::Function1Types::None

Description
    Templated placeholder function that returns an error message if called.

    This is principally useful for interfaces that expect a Function1
    but where it is not necessarily used by a particular submodel.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_None_H
#define Foam_Function1Types_None_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                            Class None Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class None
:
    public Function1<Type>
{
    // Private Member Data

        //- Context (eg, dictionary name) for the function
        string context_;


public:

    //- Runtime type information
    TypeName("none");


    // Generated Methods

        //- Default copy construct
        None(const None<Type>&) = default;

        //- No copy assignment
        void operator=(const None<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        None
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
    virtual ~None() = default;


    // Member Functions

        //- Value is independent of x
        virtual inline bool constant() const { return true; }

        //- Function cannot be evaluated
        virtual inline bool good() const { return false; }

        //- Placeholder: generates an error if called
        virtual Type value(const scalar x) const;

        //- Placeholder: generates an error if called
        virtual Type integral(const scalar x1, const scalar x2) const;

        //- Placeholder: generates an error if called
        virtual tmp<Field<Type>> value(const scalarField& x) const;

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

#ifdef NoRepository
    #include "NoneFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
