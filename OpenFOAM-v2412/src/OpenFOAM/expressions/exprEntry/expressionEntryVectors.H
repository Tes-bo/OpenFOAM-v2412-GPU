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
    Creates an expression string from common VectorSpace types

SourceFiles
    expressionEntryVectors.C

\*---------------------------------------------------------------------------*/

#ifndef exprTools_expressionEntryVectors_H
#define exprTools_expressionEntryVectors_H

#include "expressionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

#undef declareExpressionEntryType
#define declareExpressionEntryType(Type)                                      \
    /*! \brief Create an expression string from a Type */                     \
    struct Type##Entry : public exprTools::expressionEntry                    \
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


declareExpressionEntryType(vector);
declareExpressionEntryType(tensor);
declareExpressionEntryType(symmTensor);
declareExpressionEntryType(sphericalTensor);

#undef declareExpressionEntryType

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace exprTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
