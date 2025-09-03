/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "dimensionSet.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dimensionSet, 1);
}

const Foam::scalar Foam::dimensionSet::smallExponent = SMALL;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline bool checkDims
(
    const char* what,
    const dimensionSet& a,
    const dimensionSet& b
)
{
    if (a != b)
    {
        FatalErrorInFunction
            << "Different dimensions for '" << what
            << "'\n     dimensions : " << a << " != " << b << nl
            << abort(FatalError);
        return false;
    }

    return true;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSet::dimensionSet()
:
    exponents_(Zero)
{}


Foam::dimensionSet::dimensionSet
(
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles,
    const scalar current,
    const scalar luminousIntensity
)
:
    exponents_()
{
    exponents_[MASS] = mass;
    exponents_[LENGTH] = length;
    exponents_[TIME] = time;
    exponents_[TEMPERATURE] = temperature;
    exponents_[MOLES] = moles;
    exponents_[CURRENT] = current;
    exponents_[LUMINOUS_INTENSITY] = luminousIntensity;
}


Foam::dimensionSet::dimensionSet(const FixedList<scalar,7>& dims)
:
    exponents_(dims)
{}


Foam::dimensionSet::dimensionSet(const dimensionSet& ds)
:
    exponents_(ds.exponents_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dimensionSet::dimensionless() const
{
    for (const scalar val : exponents_)
    {
        // ie, mag(val) > smallExponent
        if ((val > smallExponent) || (val < -smallExponent))
        {
            return false;
        }
    }

    return true;
}


const Foam::FixedList<Foam::scalar,7>&
Foam::dimensionSet::values() const noexcept
{
    return exponents_;
}


Foam::FixedList<Foam::scalar,7>&
Foam::dimensionSet::values() noexcept
{
    return exponents_;
}


void Foam::dimensionSet::clear()
{
    exponents_ = Zero;
}


void Foam::dimensionSet::reset(const dimensionSet& ds)
{
    exponents_ = ds.exponents_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::scalar Foam::dimensionSet::operator[](const dimensionType type) const
{
    return exponents_[type];
}


Foam::scalar& Foam::dimensionSet::operator[](const dimensionType type)
{
    return exponents_[type];
}


Foam::scalar Foam::dimensionSet::operator[](const label type) const
{
    return exponents_[type];
}


Foam::scalar& Foam::dimensionSet::operator[](const label type)
{
    return exponents_[type];
}


bool Foam::dimensionSet::operator==(const dimensionSet& ds) const
{
    for (int d=0; d<nDimensions; ++d)
    {
        if
        (
            mag(exponents_[d] - ds.exponents_[d])
          > smallExponent
        )
        {
            return false;
        }
    }

    return true;
}


bool Foam::dimensionSet::operator!=(const dimensionSet& ds) const
{
    return !(operator==(ds));
}


bool Foam::dimensionSet::operator=(const dimensionSet& ds) const
{
    if (dimensionSet::checking())
    {
        checkDims("(a = b)", *this, ds);
    }

    return true;
}


bool Foam::dimensionSet::operator+=(const dimensionSet& ds) const
{
    if (dimensionSet::checking())
    {
        checkDims("(a += b)", *this, ds);
    }

    return true;
}


bool Foam::dimensionSet::operator-=(const dimensionSet& ds) const
{
    if (dimensionSet::checking())
    {
        checkDims("(a -= b)", *this, ds);
    }

    return true;
}


bool Foam::dimensionSet::operator*=(const dimensionSet& ds)
{
    reset((*this)*ds);

    return true;
}


bool Foam::dimensionSet::operator/=(const dimensionSet& ds)
{
    reset((*this)/ds);

    return true;
}


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::dimensionSet Foam::min(const dimensionSet& a, const dimensionSet& b)
{
    if (dimensionSet::checking())
    {
        checkDims("min(a, b)", a, b);
    }

    return a;
}


Foam::dimensionSet Foam::max(const dimensionSet& a, const dimensionSet& b)
{
    if (dimensionSet::checking())
    {
        checkDims("max(a, b)", a, b);
    }

    return a;
}


Foam::dimensionSet Foam::clamp(const dimensionSet& a, const dimensionSet& range)
{
    // In may cases the min/max range will be created from raw values
    // (no dimension) so accept those without error
    if (dimensionSet::checking() && !range.dimensionless())
    {
        checkDims("clamp(a, b)", a, range);
    }

    return a;
}

Foam::dimensionSet Foam::lerp(const dimensionSet& a, const dimensionSet& b)
{
    if (dimensionSet::checking())
    {
        checkDims("lerp(a, b)", a, b);
    }

    return a;
}


Foam::dimensionSet Foam::cmptMultiply
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::cmptDivide
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1/ds2;
}


Foam::dimensionSet Foam::pow(const dimensionSet& ds, const scalar p)
{
    return dimensionSet
    (
        ds[dimensionSet::MASS]*p,
        ds[dimensionSet::LENGTH]*p,
        ds[dimensionSet::TIME]*p,
        ds[dimensionSet::TEMPERATURE]*p,
        ds[dimensionSet::MOLES]*p,
        ds[dimensionSet::CURRENT]*p,
        ds[dimensionSet::LUMINOUS_INTENSITY]*p
    );
}


Foam::dimensionSet Foam::pow
(
    const dimensionSet& ds,
    const dimensionedScalar& p
)
{
    if (dimensionSet::checking() && !p.dimensions().dimensionless())
    {
        FatalErrorInFunction
            << "Exponent of pow is not dimensionless" << endl
            << abort(FatalError);
    }

    return pow(ds, p.value());
}


Foam::dimensionSet Foam::sqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


Foam::dimensionSet Foam::pow2(const dimensionSet& ds)
{
    return pow(ds, 2);
}


Foam::dimensionSet Foam::pow3(const dimensionSet& ds)
{
    return pow(ds, 3);
}


Foam::dimensionSet Foam::pow4(const dimensionSet& ds)
{
    return pow(ds, 4);
}


Foam::dimensionSet Foam::pow5(const dimensionSet& ds)
{
    return pow(ds, 5);
}


Foam::dimensionSet Foam::pow6(const dimensionSet& ds)
{
    return pow(ds, 6);
}


Foam::dimensionSet Foam::pow025(const dimensionSet& ds)
{
    return pow(ds, 0.25);
}


Foam::dimensionSet Foam::sqrt(const dimensionSet& ds)
{
    return pow(ds, 0.5);
}


Foam::dimensionSet Foam::cbrt(const dimensionSet& ds)
{
    return pow(ds, 1.0/3.0);
}


Foam::dimensionSet Foam::magSqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


Foam::dimensionSet Foam::mag(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::sign(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::pos(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::pos0(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::neg(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::neg0(const dimensionSet&)
{
    return dimless;
}


Foam::dimensionSet Foam::posPart(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::negPart(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::inv(const dimensionSet& ds)
{
    return dimensionSet
    (
        // Avoid signed zero values by subtracting from zero (with rounding)
        // instead of using a simple '-value'
        0.0-ds[dimensionSet::MASS],
        0.0-ds[dimensionSet::LENGTH],
        0.0-ds[dimensionSet::TIME],
        0.0-ds[dimensionSet::TEMPERATURE],
        0.0-ds[dimensionSet::MOLES],
        0.0-ds[dimensionSet::CURRENT],
        0.0-ds[dimensionSet::LUMINOUS_INTENSITY]
    );
}


Foam::dimensionSet Foam::trans(const dimensionSet& ds)
{
    if (dimensionSet::checking() && !ds.dimensionless())
    {
        FatalErrorInFunction
            << "Argument of trancendental function not dimensionless" << nl
            << abort(FatalError);
    }

    return ds;
}


Foam::dimensionSet Foam::atan2(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::checking())
    {
        checkDims("atan2(a, b)", ds1, ds2);
    }

    return dimless;
}


Foam::dimensionSet Foam::hypot(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::checking())
    {
        checkDims("hypot(a, b)", ds1, ds2);
    }

    return ds1;
}


Foam::dimensionSet Foam::stabilise
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    if (dimensionSet::checking())
    {
        checkDims("stabilise(a, b)", ds1, ds2);
    }

    return ds1;
}


Foam::dimensionSet Foam::transform(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::invTransform(const dimensionSet& ds)
{
    return ds;
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

Foam::dimensionSet Foam::operator~(const dimensionSet& ds)
{
    return inv(ds);
}


Foam::dimensionSet Foam::operator-(const dimensionSet& ds)
{
    return ds;
}


Foam::dimensionSet Foam::operator+
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    if (dimensionSet::checking())
    {
        checkDims("(a + b)", ds1, ds2);
    }

    return ds1;
}


Foam::dimensionSet Foam::operator-
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    if (dimensionSet::checking())
    {
        checkDims("(a - b)", ds1, ds2);
    }

    return ds1;
}


Foam::dimensionSet Foam::operator*
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet result(ds1);

    auto rhs = ds2.values().begin();

    for (scalar& val : result.values())
    {
        val += *rhs;
        ++rhs;
    }

    return result;
}


Foam::dimensionSet Foam::operator/
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet result(ds1);

    auto rhs = ds2.values().begin();

    for (scalar& val : result.values())
    {
        val -= *rhs;
        ++rhs;
    }

    return result;
}


Foam::dimensionSet Foam::operator&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::operator^
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


Foam::dimensionSet Foam::operator&&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


// ************************************************************************* //
