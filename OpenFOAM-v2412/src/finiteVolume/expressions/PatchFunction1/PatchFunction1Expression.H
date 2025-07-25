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
    Foam::PatchFunction1Types::PatchExprField

Description
    PatchFunction1 with values supplied by a parsed expression.

Usage
    Example:
    \verbatim
    <patchName>
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
                scalar(arg() > start && arg() < stop) * vector(1, 0, 0)
            #};
        }
    }
    \endverbatim

See also
    Foam::exprFixedValueFvPatchField

SourceFiles
    PatchFunction1Expression.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_expression_H
#define PatchFunction1Types_expression_H

#include "PatchFunction1.H"
#include "patchExprDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                       Class PatchExprField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PatchExprField
:
    public PatchFunction1<Type>
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        //- The expression
        expressions::exprString valueExpr_;

        //- The expression driver
        mutable expressions::patchExpr::parseDriver driver_;


public:

    //- Runtime type information
    TypeName("expression");


    // Generated Methods

        //- No copy assignment
        void operator=(const PatchExprField<Type>&) = delete;


    // Constructors

        //- Construct from patch, entry name and dictionary
        //  The patch must correspond to an fvPatch!
        PatchExprField
        (
            const polyPatch& pp,
            const word& redirectType,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct setting patch
        explicit PatchExprField
        (
            const PatchExprField<Type>& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        explicit PatchExprField(const PatchExprField<Type>& rhs);

        //- Return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return PatchFunction1<Type>::Clone(*this);
        }

        //- Return a clone, setting the patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return PatchFunction1<Type>::Clone(*this, pp);
        }


    //- Destructor
    virtual ~PatchExprField() = default;


    // Member Functions

    // Evaluation

        //- Return value.
        // The parameter 'x' is accessible as 'arg' in the function
        virtual tmp<Field<Type>> value(const scalar x) const;

        //- Is value constant (i.e. independent of x)
        virtual inline bool constant() const
        {
            return false;
        }

        //- Is value uniform (i.e. independent of coordinate)
        virtual inline bool uniform() const
        {
            return false;
        }

        //- Integrate between two values
        virtual tmp<Field<Type>> integrate
        (
            const scalar x1,
            const scalar x2
        ) const;


    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const FieldMapper& mapper);

        //- Reverse map the given PatchFunction1 onto this PatchFunction1
        virtual void rmap
        (
            const PatchFunction1<Type>& pf1,
            const labelList& addr
        );


    // IO

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFunction1Expression.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
