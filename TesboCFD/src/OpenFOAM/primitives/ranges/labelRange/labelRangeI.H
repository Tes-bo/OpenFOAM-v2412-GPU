/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::labelRange::labelRange() noexcept
:
    IntRange<label>()
{}


inline constexpr Foam::labelRange::labelRange(const label len) noexcept
:
    IntRange<label>(len)
{}


inline Foam::labelRange::labelRange
(
    const label beg,
    const label len
) noexcept
:
    IntRange<label>(beg, len)
{
    clampSize();
}


inline Foam::labelRange::labelRange
(
    const label beg,
    const label len,
    const bool adjustStart
) noexcept
:
    IntRange<label>(beg, len)
{
    if (adjustStart)
    {
        adjust(); // Eliminate negative start, adjust size accordingly
    }
    else
    {
        clampSize();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::labelRange::reset
(
    const label beg,
    const label len,
    const bool adjustStart
) noexcept
{
    reset(beg, len);

    if (adjustStart)
    {
        adjust(); // Eliminate negative start, adjust size accordingly
    }
    else
    {
        clampSize();
    }
}


// ************************************************************************* //
