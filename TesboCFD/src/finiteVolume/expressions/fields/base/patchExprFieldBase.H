/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::expressions::patchExprFieldBase

Description
    Base class for managing patches with expressions.
    The expected input supports value, gradient and mixed conditions.

Usage
    \table
        Property     | Description                          | Required | Default
        valueExpr    | expression for uniformValue          | partly  | 0
        gradientExpr | expression for uniformGradient       | partly  | 0
        fractionExpr | expression for valueFraction         | partly  | depends
    \endtable

SourceFiles
    patchExprFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_patchExprFieldBase_H
#define expressions_patchExprFieldBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "exprString.H"

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                     Class patchExprFieldBase Declaration
\*---------------------------------------------------------------------------*/

class patchExprFieldBase
{
protected:

    // Protected Types

        //- Enumeration of expected expressions
        enum expectedTypes
        {
            VALUE_TYPE = 1,
            GRADIENT_TYPE = 2,
            MIXED_TYPE = 3
        };


    // Protected Data

        //- Add debugging
        bool debug_;

        //- Slightly dodgy concept here
        bool evalOnConstruct_;

        // The expressions
        expressions::exprString valueExpr_;
        expressions::exprString gradExpr_;
        expressions::exprString fracExpr_;


private:

    // Private Member Functions

        //- Read expressions from dictionary
        void readExpressions
        (
            const dictionary& dict,
            enum expectedTypes expectedType,
            bool wantPointData = false
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Default construct
        patchExprFieldBase();

        //- Construct from dictionary
        explicit patchExprFieldBase
        (
            const dictionary& dict,
            enum expectedTypes expectedType = expectedTypes::VALUE_TYPE,
            bool wantPointData = false
        );


    // Member Functions

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
