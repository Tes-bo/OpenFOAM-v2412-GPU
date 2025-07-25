/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

inline constexpr Foam::scalarRange::scalarRange
(
    const scalarRange::rangeTypes type,
    const scalar minVal,
    const scalar maxVal
) noexcept
:
    min_(minVal),
    max_(maxVal),
    type_(type)
{}


inline constexpr Foam::scalarRange::scalarRange() noexcept
:
    scalarRange(rangeTypes::NONE, GREAT, -GREAT)
{}


inline constexpr Foam::scalarRange::scalarRange(const scalar val) noexcept
:
    scalarRange(rangeTypes::EQ, val, val)
{}


inline Foam::scalarRange::scalarRange
(
    const scalar minVal,
    const scalar maxVal
) noexcept
:
    scalarRange(rangeTypes::GE_LE, minVal, maxVal)
{
    if (maxVal < minVal)
    {
        reset();  // Inverted - explicitly mark as such
    }
    else if (equal(minVal, maxVal))
    {
        type_ = rangeTypes::EQ;
    }
}


inline constexpr Foam::scalarRange
Foam::scalarRange::ge(const scalar minVal) noexcept
{
    return scalarRange(rangeTypes::GE, minVal, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::gt(const scalar minVal) noexcept
{
    return scalarRange(rangeTypes::GT, minVal, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::ge0() noexcept
{
    return scalarRange(rangeTypes::GE, 0, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::gt0() noexcept
{
    return scalarRange(rangeTypes::GT, 0, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::le(const scalar maxVal) noexcept
{
    return scalarRange(rangeTypes::LE, -GREAT, maxVal);
}

inline constexpr Foam::scalarRange
Foam::scalarRange::lt(const scalar maxVal) noexcept
{
    return scalarRange(rangeTypes::LT, -GREAT, maxVal);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::zero_one() noexcept
{
    return scalarRange(rangeTypes::GE_LE, 0, 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::scalarRange::reset() noexcept
{
    min_ = GREAT;
    max_ = -GREAT;
    type_ = rangeTypes::NONE;
}


inline bool Foam::scalarRange::empty() const noexcept
{
    return (type_ == rangeTypes::NONE);
}


inline bool Foam::scalarRange::good() const noexcept
{
    return (type_ != rangeTypes::NONE);
}


inline bool Foam::scalarRange::single() const noexcept
{
    return (type_ == rangeTypes::EQ);
}


inline Foam::scalar Foam::scalarRange::value() const
{
    switch (type_)
    {
        case rangeTypes::EQ: // For equals, min and max are identical
        case rangeTypes::GE:
        case rangeTypes::GT:
            return min_;

        case rangeTypes::LE:
        case rangeTypes::LT:
            return max_;

        case rangeTypes::GE_LE:
            // Multiply before adding to avoid possible overflow
            return (0.5 * min_) + (0.5 * max_);

        default:  // NONE or ALWAYS. Zero is reasonable dummy value.
            return 0;
    }
}


inline bool Foam::scalarRange::contains(const scalar val) const
{
    switch (type_)
    {
        case rangeTypes::EQ:    return equal(val, min_);
        case rangeTypes::GE:    return (val >= min_);
        case rangeTypes::GT:    return (val > min_);
        case rangeTypes::LE:    return (val <= max_);
        case rangeTypes::LT:    return (val < max_);
        case rangeTypes::GE_LE: return (val >= min_ && val <= max_);
        case rangeTypes::ALWAYS: return true;
        default:    return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline constexpr bool
Foam::scalarRange::operator==(const scalarRange& rhs) const noexcept
{
    return (type_ == rhs.type_ && min_ == rhs.min_ && max_ == rhs.max_);
}


inline constexpr bool
Foam::scalarRange::operator!=(const scalarRange& rhs) const noexcept
{
    return !(*this == rhs);
}


// ************************************************************************* //
