/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "uint64.H"
#include "parsing.H"
#include "IOstreams.H"
#include <cinttypes>
#include <cmath>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint64_t Foam::readUint64(const char* buf)
{
    char *endptr = nullptr;
    errno = 0;
    const uintmax_t parsed = ::strtoumax(buf, &endptr, 10);

    const uint64_t val = uint64_t(parsed);

    const parsing::errorType err =
    (
        (parsed > UINT64_MAX)
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


bool Foam::readUint64(const char* buf, uint64_t& val)
{
    char *endptr = nullptr;
    errno = 0;
    const uintmax_t parsed = ::strtoumax(buf, &endptr, 10);

    val = uint64_t(parsed);

    return
    (
        (parsed > UINT64_MAX)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


uint64_t Foam::readUint64(Istream& is)
{
    uint64_t val(0);
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, uint64_t& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get uint64"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        val = uint64_t(t.labelToken());
    }
    else if (t.isScalar())
    {
        const scalar sval(t.scalarToken());
        const uintmax_t parsed = uintmax_t(std::round(sval));
        val = 0 + uint64_t(parsed);

        // Accept integral floating-point values.
        // Eg, from string expression evaluation (#1696)

        if ((sval < -1e-4) || parsed > UINT64_MAX)
        {
            FatalIOErrorInFunction(is)
                << "Expected label (uint64), value out-of-range "
                << t.info()
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
        else if (1e-4 < std::abs(sval - scalar(parsed)))
        {
            FatalIOErrorInFunction(is)
                << "Expected label (uint64), found non-integral value "
                << t.info()
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (uint64), found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const uint64_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


#ifdef __APPLE__
Foam::Ostream& Foam::operator<<(Ostream& os, const unsigned long val)
{
    os << uint64_t(val);
    return os;
}
#endif


// ************************************************************************* //
