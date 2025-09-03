/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 Bernhard Gschaider
    Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Creates an expression string from various standard dimensioned types.

SourceFiles
    expressionEntryDimensioned.C

\*---------------------------------------------------------------------------*/

#ifndef exprTools_expressionEntryDimensioned_H
#define exprTools_expressionEntryDimensioned_H

#include "expressionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

#undef declareExpressionEntryType
#define declareExpressionEntryType(DimType)                                   \
    /*! \brief Create an expression string from a DimType as its value */     \
    struct DimType##Entry : public exprTools::expressionEntry                 \
    {                                                                         \
        /*! \brief Entry to string */                                         \
        static string evaluate(const entry& e);                               \
                                                                              \
        /*! \brief Entry to string */                                         \
        virtual string toExpr(const entry& e) const                           \
        {                                                                     \
            return evaluate(e);                                               \
        }                                                                     \
    };


declareExpressionEntryType(dimensionedScalar);
declareExpressionEntryType(dimensionedVector);
declareExpressionEntryType(dimensionedTensor);
declareExpressionEntryType(dimensionedSymmTensor);
declareExpressionEntryType(dimensionedSphericalTensor);

#undef declareExpressionEntryType

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace exprTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
