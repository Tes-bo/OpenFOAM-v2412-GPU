/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::find
(
    const T& val,
    label pos,
    label len
) const
{
    if (pos >= 0 && pos < addr_.size())
    {
        // Change sub-length to (one-past) end position
        // len == -1 (like std::string::npos) - search until end

        if (len > 0) len += pos;
        if (len < 0 || len > addr_.size())
        {
            len = addr_.size();
        }

        const T* const vals = values_.begin();

        while (pos < len)
        {
            if (vals[addr_[pos]] == val)
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::rfind
(
    const T& val,
    label pos
) const
{
    // pos == -1 (like std::string::npos) - search from end

    if (pos < 0 || pos >= addr_.size())
    {
        pos = addr_.size()-1;
    }

    const T* const vals = values_.begin();

    while (pos >= 0)
    {
        if (vals[addr_[pos]] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


// ************************************************************************* //
