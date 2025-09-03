/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "expressionEntryVectors.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    vectorEntry,
    empty,
    vector
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    tensorEntry,
    empty,
    tensor
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    symmTensorEntry,
    empty,
    symmTensor
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    sphericalTensorEntry,
    empty,
    sphericalTensor
);

} // End namespace exprTools
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#undef defineExpressionEntryType
#define defineExpressionEntryType(Type)                                       \
    Foam::string Foam::exprTools::Type##Entry::evaluate(const entry& e)       \
    {                                                                         \
        return toExprStr<Type>(e.stream());                                   \
    }


defineExpressionEntryType(vector);
defineExpressionEntryType(tensor);
defineExpressionEntryType(symmTensor);
defineExpressionEntryType(sphericalTensor);

#undef defineExpressionEntryType

// ************************************************************************* //
