/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "uint8.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<uint8_t>::typeName = "uint8";
const char* const Foam::pTraits<uint8_t>::componentNames[] = { "" };

const uint8_t Foam::pTraits<uint8_t>::zero = 0;
const uint8_t Foam::pTraits<uint8_t>::one = 1;
const uint8_t Foam::pTraits<uint8_t>::min = 0;
const uint8_t Foam::pTraits<uint8_t>::max = UINT8_MAX;
const uint8_t Foam::pTraits<uint8_t>::rootMin = 0;
const uint8_t Foam::pTraits<uint8_t>::rootMax = UINT8_MAX;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pTraits<uint8_t>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint8_t Foam::readUint8(Istream& is)
{
    uint8_t val(0);
    is >> val;

    return val;
}

// IO operators are identical to direction, which is uint8_t


// ************************************************************************* //
