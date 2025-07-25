/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "DynamicList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::label Foam::DynamicList<T, SizeMin>::removeElements
(
    const labelRange& slice
)
{
    if (!slice.good())
    {
        // No-op
        return 0;
    }

    // Note: already checked for valid begin_value before
    if (slice.end_value() >= this->size())
    {
        // Remove entire tail
        this->resize(slice.begin_value());
    }
    else
    {
        // Move all trailing elements down into space previously
        // occupied by the slice. Truncate after
        std::move
        (
            this->begin(slice.end_value()),
            this->end(),
            this->begin(slice.begin_value())
        );

        this->resize(this->size() - slice.size());
    }

    return slice.size();
}


template<class T, int SizeMin>
Foam::label Foam::DynamicList<T, SizeMin>::subsetElements
(
    const labelRange& slice
)
{
    if (slice.begin_value() > 0)
    {
        // Move elements down.
        // Since begin_value > 0, the initial destination is non-overlapping
        std::move
        (
            this->begin(slice.begin_value()),
            this->begin(slice.end_value()),
            this->begin()
        );
    }

    // Don't need min size, since slice size was already checked before
    resize(slice.size());
    return this->size();
}


// ************************************************************************* //
