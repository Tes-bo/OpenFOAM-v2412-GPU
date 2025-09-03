/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "HashOps.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelHashSet Foam::HashSetOps::used(const bitSet& select)
{
    labelHashSet output;

    if (select.any())
    {
        output.reserve(select.count());

        for (label i = select.find_first(); i >= 0; i = select.find_next(i))
        {
            output.insert(i);
        }
    }

    return output;
}


Foam::labelHashSet Foam::HashSetOps::used(const UList<bool>& select)
{
    const label len = select.size();

    // No idea of the sparseness, just assume 1/8
    labelHashSet output(len/4);

    for (label i = 0; i < len; ++i)
    {
        if (select[i])
        {
            output.insert(i);
        }
    }

    return output;
}


Foam::bitSet Foam::HashSetOps::bitset(const labelHashSet& locations)
{
    bitSet output;
    output.setMany(locations.begin(), locations.end());

    return output;
}


Foam::List<bool> Foam::HashSetOps::bools(const labelHashSet& locations)
{
    auto const max = std::max_element(locations.begin(), locations.end());
    const label len = (max != locations.end() ? (1 + *max) : 0);

    if (len <= 0)
    {
        return List<bool>();
    }

    List<bool> output(len, false);

    for (const label i : locations)
    {
        if (i >= 0)
        {
            output[i] = true;
        }
    }

    return output;
}


// ************************************************************************* //
