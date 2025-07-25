/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::Function1Types::Function1Expression

Description
    Function1 with values supplied by a parsed expression.

Usage
    Example:
    \verbatim
    \<patchName\>
    {
        type            uniformFixedValue;
        uniformValue
        {
            type            expression;

            // optional variables for use within the expression
            variables
            (
                "start = 0.5"
                "stop  = 1"
            );

            // A step function
            expression
            #{
                mag(arg() > start && arg() < stop) * vector(1, 0, 0)
            #};
        }
    }
    \endverbatim

Very much like '\#eval' but runtime

See also
    Foam::exprFixedValueFvPatchField

SourceFiles
    Function1Expression.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_expression_H
#define Foam_Function1Types_expression_H

#include "Function1.H"
#include "fieldExprDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                     Class Function1Expression Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1Expression
:
    public Function1<Type>
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        //- The expression
        expressions::exprString valueExpr_;

        //- The expression driver
        mutable expressions::fieldExprDriver driver_;


public:

    //- Runtime type information
    TypeName("expression");


    // Generated Methods

        //- No copy assignment
        void operator=(const Function1Expression<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Function1Expression
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit Function1Expression(const Function1Expression<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Function1Expression() = default;


    // Member Functions

        //- Return value.
        //  The parameter 'x' is accessible as 'arg' in the function
        virtual Type value(const scalar x) const;

        //- Integrate between two values
        virtual Type integrate
        (
            const scalar x1,
            const scalar x2
        ) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1Expression.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
