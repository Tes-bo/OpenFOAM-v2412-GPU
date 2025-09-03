/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

inline Foam::boolVector::boolVector()
:
    FixedList<bool, 3>(false)
{}


inline Foam::boolVector::boolVector(const bool val)
:
    FixedList<bool, 3>(val)
{}


inline Foam::boolVector::boolVector
(
    const bool vx,
    const bool vy,
    const bool vz
)
{
    x() = vx;
    y() = vy;
    z() = vz;
}


inline Foam::boolVector::boolVector(Istream& is)
:
    FixedList<bool, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::boolVector::all() const noexcept
{
    return (x() && y() && z());
}


inline bool Foam::boolVector::any() const noexcept
{
    return (x() || y() || z());
}


inline bool Foam::boolVector::none() const noexcept
{
    return !any();
}


inline unsigned int Foam::boolVector::count(const bool on) const
{
    unsigned int total = 0;

    for (const bool val : *this)
    {
        if (val) ++total;
    }

    if (!on)
    {
        // Return the number of bits that are OFF.
        return (3u - total);
    }

    return total;
}


inline void Foam::boolVector::flip()
{
    x() = !x();
    y() = !y();
    z() = !z();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::boolVector::operator=(const bool value)
{
    FixedList<bool, 3>::operator=(value);
}


// ************************************************************************* //
