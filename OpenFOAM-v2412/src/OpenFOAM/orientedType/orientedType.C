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

#include "orientedType.H"
#include "dictionary.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::orientedType::orientedOption
>
Foam::orientedType::orientedOptionNames
({
    { orientedOption::UNKNOWN, "unknown" },
    { orientedOption::ORIENTED, "oriented" },
    { orientedOption::UNORIENTED, "unoriented" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::orientedType::checkType
(
    const orientedType& a,
    const orientedType& b
) noexcept
{
    return
    (
        (a.oriented() == b.oriented())
     || (a.oriented() == orientedOption::UNKNOWN)
     || (b.oriented() == orientedOption::UNKNOWN)
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline bool checkTypes
(
    const char* what,
    const orientedType& a,
    const orientedType& b
)
{
    // ie, checkType(a,b)
    const bool ok
    (
        (a.oriented() == b.oriented())
     || (a.oriented() == orientedType::orientedOption::UNKNOWN)
     || (b.oriented() == orientedType::orientedOption::UNKNOWN)
    );

    if (!ok)
    {
        FatalErrorInFunction
            << what << " : undefined for "
            << orientedType::orientedOptionNames[a.oriented()] << " and "
            << orientedType::orientedOptionNames[b.oriented()] << " types"
            << abort(FatalError);
    }

    return ok;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orientedType::orientedType(Istream& is)
:
    oriented_(orientedOptionNames.read(is))
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::orientedType::read(const dictionary& dict)
{
    oriented_ = orientedOptionNames.getOrDefault
    (
        "oriented",
        dict,
        orientedOption::UNKNOWN,
        true  // Failsafe behaviour
    );
}


bool Foam::orientedType::writeEntry(Ostream& os) const
{
    const bool output = (oriented_ == ORIENTED);

    if (output)
    {
        os.writeEntry("oriented", orientedOptionNames[oriented_]);
    }

    return output;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::orientedType::operator+=(const orientedType& ot)
{
    // Set the oriented status if it was unknown
    if (oriented_ == UNKNOWN)
    {
        oriented_ = ot.oriented();
    }

    (void) checkTypes("Operator +=", *this, ot);
}


void Foam::orientedType::operator-=(const orientedType& ot)
{
    // Set the oriented status if it was unknown
    if (oriented_ == orientedOption::UNKNOWN)
    {
        oriented_ = ot.oriented();
    }

    (void) checkTypes("Operator -=", *this, ot);
}


void Foam::orientedType::operator*=(const orientedType& ot)
{
    setOriented(is_oriented() != ot.is_oriented());
}


void Foam::orientedType::operator/=(const orientedType& ot)
{
    setOriented(is_oriented() != ot.is_oriented());
}


void Foam::orientedType::operator*=(const scalar s)
{
    // No change to oriented_ flag
}


void Foam::orientedType::operator/=(const scalar s)
{
    // No change to oriented_ flag
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::orientedType Foam::min(const orientedType& a, const orientedType& b)
{
    (void) checkTypes("Function min", a, b);
    return a;
}


Foam::orientedType Foam::max(const orientedType& a, const orientedType& b)
{
    (void) checkTypes("Function max", a, b);
    return a;
}


Foam::orientedType Foam::lerp(const orientedType& a, const orientedType& b)
{
    (void) checkTypes("Function lerp", a, b);

    // Behaves like addition
    return orientedType(a.is_oriented() || b.is_oriented());
}


Foam::orientedType Foam::cmptMultiply
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::cmptDivide
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::cmptAv(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow(const orientedType& ot, const scalar p)
{
    // Undefined???
    // - only defined for integers where:
    //   - odd powers = oriented_ = yes (if ot is oriented)
    //   - even powers = oriented_ = no
    return ot;
}


Foam::orientedType Foam::sqr(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow2(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow3(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow4(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow5(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow6(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow025(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::sqrt(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::cbrt(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::magSqr(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType  Foam::mag(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType  Foam::sign(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::pos(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::pos0(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::neg(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::neg0(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::posPart(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::negPart(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::inv(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::trans(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::atan2
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    (void) checkTypes("Function atan2", ot1, ot2);

    return ot1;
}


Foam::orientedType Foam::hypot
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    (void) checkTypes("Function hypot", ot1, ot2);

    return ot1;
}


Foam::orientedType Foam::transform(const orientedType& ot)
{
    return ot;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, orientedType& ot)
{
    ot.oriented_ = orientedType::orientedOptionNames.read(is);

    is.check(FUNCTION_NAME);

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const orientedType& ot)
{
    os << orientedType::orientedOptionNames[ot.oriented()];

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::orientedType Foam::operator+
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    (void) checkTypes("Operator +", ot1, ot2);

    return orientedType(ot1.is_oriented() || ot2.is_oriented());
}


Foam::orientedType Foam::operator-(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::operator-
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    (void) checkTypes("Operator -", ot1, ot2);

    return orientedType(ot1.is_oriented() || ot2.is_oriented());
}


Foam::orientedType Foam::operator*(const scalar s, const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::operator/(const orientedType& ot, const scalar s)
{
    return ot;
}


Foam::orientedType Foam::operator/
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator*
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator^
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return orientedType(ot1.is_oriented() != ot2.is_oriented());
}


Foam::orientedType Foam::operator&
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator&&
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return orientedType(false);
}


// ************************************************************************* //
