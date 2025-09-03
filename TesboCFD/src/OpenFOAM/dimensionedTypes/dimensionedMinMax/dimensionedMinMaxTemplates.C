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

#include "dimensionedMinMax.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


template<class T>
Foam::dimensioned<Foam::MinMax<T>> Foam::makeDimensionedMinMax
(
    const word& name,
    const dimensionSet& dims,
    const MinMax<T>& values,
    const dictionary& dict,
    const word& minName,
    const word& maxName
)
{
    // Normal construction with optional entry

    dimensioned<MinMax<T>> range(name, dims, values, dict);

    // Optional min specification
    if (!minName.empty())
    {
        dimensioned<T> minVal(minName, dims, values.min(), dict);
        range.dimensions() += minVal.dimensions();
        range.value().min() = minVal.value();
    }

    // Optional max specification
    if (!maxName.empty())
    {
        dimensioned<T> maxVal(maxName, dims, values.max(), dict);
        range.dimensions() += maxVal.dimensions();
        range.value().max() = maxVal.value();
    }

    return range;
}


// ************************************************************************* //
