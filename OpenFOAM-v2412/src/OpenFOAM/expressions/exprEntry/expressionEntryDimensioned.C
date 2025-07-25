/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "expressionEntryDimensioned.H"
#include "primitiveEntry.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "dimensionedTensor.H"
#include "dimensionedSymmTensor.H"
#include "dimensionedSphericalTensor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    dimensionedScalarEntry,
    empty,
    dimensionedScalar
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    dimensionedVectorEntry,
    empty,
    dimensionedVector
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    dimensionedTensorEntry,
    empty,
    dimensionedTensor
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    dimensionedSymmTensorEntry,
    empty,
    dimensionedSymmTensor
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    dimensionedSphericalTensorEntry,
    empty,
    dimensionedSphericalTensor
);

} // End namespace exprTools
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#undef defineExpressionEntryType
#define defineExpressionEntryType(DimType)                                    \
    Foam::string Foam::exprTools::DimType##Entry::evaluate(const entry& e)    \
    {                                                                         \
        DimType dt(dynamicCast<const primitiveEntry>(e));                     \
        return toExprStr<DimType::value_type>(dt.value());                    \
    }


Foam::string Foam::exprTools::dimensionedScalarEntry::evaluate(const entry& e)
{
    dimensionedScalar dt(dynamicCast<const primitiveEntry>(e));
    return std::to_string(dt.value());
}


defineExpressionEntryType(dimensionedVector);
defineExpressionEntryType(dimensionedTensor);
defineExpressionEntryType(dimensionedSymmTensor);
defineExpressionEntryType(dimensionedSphericalTensor);

#undef defineExpressionEntryType

// ************************************************************************* //
