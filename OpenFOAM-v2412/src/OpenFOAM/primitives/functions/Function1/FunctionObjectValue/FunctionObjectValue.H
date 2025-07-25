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
    Foam::Function1Types::FunctionObjectValue

Description
    Returns a value retrieved from a function object result.

    Usage:
    \verbatim
        <entryName> functionObjectValue;
        <entryName>Coeffs
        {
            functionObject          <name>;
            functionObjectResult    <function object result field>;
            defaultValue            <value>; // optional
        }
    \endverbatim

Note
- does not implement integrate()

SourceFiles
    FunctionObjectValue.C
    FunctionObjectValueI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_FunctionObjectValue_H
#define Foam_Function1Types_FunctionObjectValue_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                     Class FunctionObjectValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FunctionObjectValue
:
    public Function1<Type>
{
    // Private Data

        //- Function object name
        word foName_;

        //- Function object result name
        word foResultName_;

        //- Default value when function object result is unavailable
        Type defaultValue_;

        //- Flag to indicate that default value is available
        bool haveDefaultValue_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    //- Runtime type information
    TypeName("functionObjectValue");


    // Generated Methods

        //- No copy assignment
        void operator=(const FunctionObjectValue<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        FunctionObjectValue
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit FunctionObjectValue(const FunctionObjectValue<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~FunctionObjectValue() = default;


    // Member Functions

        //- Return value at current time
        virtual inline Type value(const scalar /*unused*/) const;

        //- Integrate between two (scalar) values. Not implemented!
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FunctionObjectValueI.H"

#ifdef NoRepository
    #include "FunctionObjectValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
