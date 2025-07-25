/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "labelList.H"
#include <numeric>

// * * * * * * * * * * * * * * Global Data Members * * * * * * * * * * * * * //

// This was deprecated (2019-02) in favour of labelList::null()
const Foam::labelList Foam::emptyLabelList;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::identity(Foam::UList<int32_t>& map, int32_t start)
{
    std::iota(map.begin(), map.end(), start);
}


void Foam::identity(Foam::UList<int64_t>& map, int64_t start)
{
    std::iota(map.begin(), map.end(), start);
}


Foam::labelList Foam::identity(const label len, label start)
{
    labelList map(len);
    std::iota(map.begin(), map.end(), start);

    return map;
}


// ************************************************************************* //
