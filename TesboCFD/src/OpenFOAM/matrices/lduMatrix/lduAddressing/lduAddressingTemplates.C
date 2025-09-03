/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "lduAddressing.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lduAddressing::map
(
    const UList<Type>& faceVals,
    List<Type>& vals
) const
{
    const auto& offsets = losortStartAddr();
    const auto& indexToFace = losortAddr();

    const label n = offsets.size()-1;

    vals.resize_nocopy(faceVals.size());

    for (label celli = 0; celli < n; celli++)
    {
        const label start = offsets[celli];
        const label end = offsets[celli+1];

        for (label i = start; i < end; i++)
        {
            const label facei = indexToFace[i];
            vals[i] = faceVals[facei];
        }
    }
}


// ************************************************************************* //
