/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "int32.H"
#include "error.H"
#include "parsing.H"
#include "IOstreams.H"
#include <cinttypes>
#include <cmath>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int32_t Foam::readInt32(const char* buf)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    const int32_t val = int32_t(parsed);

    const parsing::errorType err =
    (
        (parsed < INT32_MIN || parsed > INT32_MAX)
      ? parsing::errorType::RANGE
      : parsing::checkConversion(buf, endptr)
    );

    if (err != parsing::errorType::NONE)
    {
        FatalIOErrorInFunction("unknown")
            << parsing::errorNames[err] << " '" << buf << "'"
            << exit(FatalIOError);
    }

    return val;
}


bool Foam::readInt32(const char* buf, int32_t& val)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    val = int32_t(parsed);

    return
    (
        (parsed < INT32_MIN || parsed > INT32_MAX)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


int32_t Foam::readInt32(Istream& is)
{
    int32_t val(0);
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, int32_t& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get int32"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    // Accept separated '-' (or '+') while expecting a number.
    // This can arise during dictionary expansions (Eg, -$value)

    char prefix = 0;
    if (t.isPunctuation())
    {
        prefix = t.pToken();
        if (prefix == token::PLUS || prefix == token::MINUS)
        {
            is >> t;
        }
    }

    if (t.isLabel())
    {
        val = int32_t
        (
            (prefix == token::MINUS)
          ? (0 - t.labelToken())
          : t.labelToken()
        );
    }
    else if (t.isScalar())
    {
        const scalar sval
        (
            (prefix == token::MINUS)
          ? (0 - t.scalarToken())
          : t.scalarToken()
        );
        const intmax_t parsed = intmax_t(std::round(sval));
        val = 0 + int32_t(parsed);

        // Accept integral floating-point values.
        // Eg, from string expression evaluation (#1696)

        if (parsed < INT32_MIN || parsed > INT32_MAX)
        {
            FatalIOErrorInFunction(is)
                << "Expected integral (int32), value out-of-range "
                << t.info()
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
        else if (1e-4 < std::abs(sval - scalar(parsed)))
        {
            FatalIOErrorInFunction(is)
                << "Expected integral (int32), found non-integral value "
                << t.info()
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (int32), found ";
        if (prefix == token::PLUS || prefix == token::MINUS)
        {
            FatalIOError << '\'' << prefix << "' followed by ";
        }
        FatalIOError << t.info() << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int32_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


#if (__SIZEOF_LONG__ == 4)
Foam::Istream& Foam::operator>>(Istream& is, long& val)
{
    return operator>>(is, reinterpret_cast<int32_t&>(val));
}

Foam::Ostream& Foam::operator<<(Ostream& os, const long val)
{
    return (os << int32_t(val));
}
#endif


// ************************************************************************* //
