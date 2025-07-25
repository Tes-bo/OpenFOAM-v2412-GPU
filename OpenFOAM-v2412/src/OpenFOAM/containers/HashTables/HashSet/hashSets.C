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

#include "HashSet.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::label Foam::min(const labelHashSet& set, label minValue)
{
    for (const label val : set)
    {
        if (minValue > val)
        {
            minValue = val;
        }
    }

    return minValue;
}


Foam::label Foam::max(const labelHashSet& set, label maxValue)
{
    for (const label val : set)
    {
        if (maxValue < val)
        {
            maxValue = val;
        }
    }

    return maxValue;
}


Foam::MinMax<Foam::label> Foam::minMax(const labelHashSet& set)
{
    MinMax<label> result;

    for (const label val : set)
    {
        result += val;
    }

    return result;
}


// ************************************************************************* //
