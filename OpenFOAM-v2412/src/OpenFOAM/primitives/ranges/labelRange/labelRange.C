/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "labelRange.H"
#include "List.H"
#include "MinMax.H"
#include "Pair.H"
#include <numeric>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    int labelRange::debug(debug::debugSwitch("labelRange", 0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelRange::labelRange(const MinMax<label>& range) noexcept
:
    labelRange()
{
    if (range.min() <= range.max())
    {
        // max is inclusive, so size with +1. Hope for no overflow
        start() = range.min();
        size()  = 1 + (range.max() - range.min());
    }
}


Foam::labelRange::labelRange(const Pair<label>& start_end) noexcept
:
    labelRange()
{
    if (start_end.first() <= start_end.second())
    {
        // second is exclusive, so size directly. Hope for no overflow
        start() = start_end.first();
        size()  = (start_end.second() - start_end.first());
    }
}


Foam::labelRange::labelRange(Istream& is)
:
    labelRange()
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::label> Foam::labelRange::labels() const
{
    List<label> result;

    if (this->size() > 0)
    {
        result.resize(this->size());
        std::iota(result.begin(), result.end(), this->start());
    }

    return result;
}


void Foam::labelRange::adjust() noexcept
{
    if (this->start() < 0)
    {
        if (this->size() > 0)
        {
            // Decrease size accordingly
            this->size() += this->start();
        }
        this->start() = 0;
    }
    clampSize();
}


bool Foam::labelRange::overlaps(const labelRange& range, bool touches) const
{
    const label extra = touches ? 1 : 0;

    return
    (
        this->size() && range.size()
     &&
        (
            (
                range.min() >= this->min()
             && range.min() <= this->max() + extra
            )
         ||
            (
                this->min() >= range.min()
             && this->min() <= range.max() + extra
            )
        )
    );
}


Foam::labelRange Foam::labelRange::join(const labelRange& range) const
{
    // Trivial cases first
    if (!this->size())
    {
        return *this;
    }
    else if (!range.size())
    {
        return range;
    }

    const label lower = Foam::min(this->min(), range.min());
    const label upper = Foam::max(this->max(), range.max());
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    labelRange newRange(lower, total);
    newRange.clampSize();

    return newRange;
}


Foam::labelRange Foam::labelRange::subset(const labelRange& range) const
{
    const label lower = Foam::max(this->min(), range.min());
    const label upper = Foam::min(this->max(), range.max());
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }

    return labelRange();
}


Foam::labelRange Foam::labelRange::subset
(
    const label start,
    const label size
) const
{
    const label lower = Foam::max(this->min(), start);
    const label upper = Foam::min(this->max(), start+Foam::max(0,size-1));
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }

    return labelRange();
}


Foam::labelRange Foam::labelRange::subset0(const label size) const
{
    const label lower = Foam::max(this->min(), 0);
    const label upper = Foam::min(this->max(), Foam::max(0,size-1));
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }

    return labelRange();
}


// ************************************************************************* //
