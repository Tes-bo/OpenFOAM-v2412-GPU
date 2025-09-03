/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Switch.H"
#include "scalar.H"
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace
{
static_assert
(
    Foam::Switch::INVALID+1 == 9,
    "Switch::switchType does not have 9 entries"
);

//- The names corresponding to the Switch::switchType enumeration.
//  Includes extra entries for "invalid".
static const char* names[9] =
{
    "false", "true",
    "no",    "yes",
    "off",   "on",
    "none",  "any",
    "invalid"  //< Output representation only
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
template<class OS>
static OS& printTokenError(OS& os, const token& tok)
{
    if (!tok.good())
    {
        os  << "Bad token - could not get bool/switch" << nl;
    }
    else if (tok.isWord())
    {
        os  << "Expected true/false, on/off... found "
            << tok.wordToken() << nl;
    }
    else
    {
        os  << "Wrong token - expected bool/switch, found "
            << tok.info() << nl;
    }

    return os;
}
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Switch::switchType Foam::Switch::parse
(
    const std::string& str,
    const bool failOnError
)
{
    switch (str.size())
    {
        case 1: // (0|1|f|t|n|y)
        {
            switch (str[0])
            {
                case '0': case 'f': return switchType::FALSE;
                case '1': case 't': return switchType::TRUE;
                case 'n': return switchType::NO;
                case 'y': return switchType::YES;
            }
            break;
        }
        case 2: // (no|on)
        {
            if (str == names[switchType::NO]) return switchType::NO;
            if (str == names[switchType::ON]) return switchType::ON;
            break;
        }
        case 3: // (off|yes|any)
        {
            if (str == names[switchType::OFF]) return switchType::OFF;
            if (str == names[switchType::YES]) return switchType::YES;
            if (str == names[switchType::ANY]) return switchType::ANY;
            break;
        }
        case 4: // (none|true)
        {
            if (str == names[switchType::NONE]) return switchType::NONE;
            if (str == names[switchType::TRUE]) return switchType::TRUE;
            break;
        }
        case 5: // (false)
        {
            if (str == names[switchType::FALSE]) return switchType::FALSE;
            break;
        }
    }

    if (failOnError)
    {
        FatalErrorInFunction
            << "Unknown switch " << str << nl
            << abort(FatalError);
    }

    return switchType::INVALID;
}


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

const char* Foam::Switch::name(const bool b) noexcept
{
    return names[(b ? 1 : 0)];
}


Foam::Switch Foam::Switch::find(const std::string& str)
{
    return Switch(parse(str, false));  // failOnError=false
}


bool Foam::Switch::contains(const std::string& str)
{
    return (switchType::INVALID != parse(str, false));  // failOnError=false
}


Foam::Switch Foam::Switch::getOrDefault
(
    const word& key,
    const dictionary& dict,
    const Switch deflt
)
{
    return dict.getOrDefault<Switch>(key, deflt, keyType::LITERAL);
}


Foam::Switch Foam::Switch::getOrAddToDict
(
    const word& key,
    dictionary& dict,
    const Switch deflt
)
{
    return dict.getOrAdd<Switch>(key, deflt, keyType::LITERAL);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Switch::Switch(const std::string& str)
:
    value_(parse(str, true))
{}


Foam::Switch::Switch(const char* str)
:
    value_(parse(str, true))
{}


Foam::Switch::Switch(const std::string& str, bool allowBad)
:
    value_(parse(str, !allowBad))
{}


Foam::Switch::Switch(const char* str, bool allowBad)
:
    value_(parse(str, !allowBad))
{}


Foam::Switch::Switch(const float val, const float tol)
:
    value_((mag(val) > tol) ? switchType::TRUE : switchType::FALSE)
{}


Foam::Switch::Switch(const double val, const double tol)
:
    value_((mag(val) > tol) ? switchType::TRUE : switchType::FALSE)
{}


Foam::Switch::Switch(const token& tok)
:
    value_(switchType::INVALID)
{
    if (tok.good())
    {
        if (tok.isBool())
        {
            (*this) = tok.boolToken();
        }
        else if (tok.isLabel())
        {
            (*this) = bool(tok.labelToken());
        }
        else if (tok.isWord())
        {
            value_ = parse(tok.wordToken(), false);  // failOnError=false
        }
    }
}


Foam::Switch::Switch
(
    const word& key,
    const dictionary& dict
)
:
    value_(switchType::INVALID)
{
    const token tok(dict.get<token>(key, keyType::LITERAL));

    Switch sw(tok);

    if (sw.good())
    {
        (*this) = sw;
    }
    else
    {
        printTokenError(FatalIOErrorInFunction(dict), tok)
            << exit(FatalIOError);
    }
}


Foam::Switch::Switch
(
    const word& key,
    const dictionary& dict,
    const Switch deflt,
    const bool warnOnly
)
:
    value_(deflt.value_)
{
    token tok;

    if (dict.readIfPresent<token>(key, tok, keyType::LITERAL))
    {
        Switch sw(tok);

        if (sw.good())
        {
            (*this) = sw;
        }
        else if (warnOnly)
        {
            printTokenError(IOWarningInFunction(dict), tok)
                << "using default " << deflt.c_str() << endl;
        }
        else
        {
            printTokenError(FatalIOErrorInFunction(dict), tok)
                << exit(FatalIOError);
        }
    }
}


Foam::Switch::Switch(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::good() const noexcept
{
    return (value_ < switchType::INVALID);
}


Foam::Switch::switchType Foam::Switch::type() const noexcept
{
    return switchType(value_);
}


void Foam::Switch::negate() noexcept
{
    if (value_ < switchType::INVALID)
    {
        // Toggle final bit. So NO <-> YES, OFF <-> ON ...
        value_ ^= 0x1;
    }
}


const char* Foam::Switch::c_str() const noexcept
{
    return names[(value_ & 0x0F)];
}


std::string Foam::Switch::str() const
{
    return names[(value_ & 0x0F)];
}


bool Foam::Switch::readIfPresent
(
    const word& key,
    const dictionary& dict
)
{
    return dict.readIfPresent<Switch>(key, *this, keyType::LITERAL);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, Switch& sw)
{
    token tok(is);

    sw = Switch(tok);

    if (sw.good())
    {
        is.check(FUNCTION_NAME);
    }
    else
    {
        printTokenError(FatalIOErrorInFunction(is), tok)
            << exit(FatalIOError);
        is.setBad();
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const Switch sw)
{
    os << sw.c_str();
    return os;
}


// ************************************************************************* //
