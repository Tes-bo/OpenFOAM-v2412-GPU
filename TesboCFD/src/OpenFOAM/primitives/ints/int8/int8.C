/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "int8.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<int8_t>::typeName = "int8";
const char* const Foam::pTraits<int8_t>::componentNames[] = { "" };

const int8_t Foam::pTraits<int8_t>::zero = 0;
const int8_t Foam::pTraits<int8_t>::one = 1;
const int8_t Foam::pTraits<int8_t>::min = INT8_MIN;
const int8_t Foam::pTraits<int8_t>::max = INT8_MAX;
const int8_t Foam::pTraits<int8_t>::rootMin = INT8_MIN;
const int8_t Foam::pTraits<int8_t>::rootMax = INT8_MAX;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pTraits<int8_t>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// int8_t Foam::readInt8(Istream& is)
// {
//     int8_t val(0);
//     is >> val;
//
//     return val;
// }


Foam::Istream& Foam::operator>>(Istream& is, int8_t& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get int8"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    // Accept separated '-' (or '+') while expecting a number.
    // This can arise during dictionary expansions (Eg, -$value)

    if (t.isLabel())
    {
        val = int8_t(t.labelToken());
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (int8), found "
            << t.info() << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int8_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
