/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "gradingDescriptors.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Prepend a value by shifting contents
template<class T>
static void prependList(List<T>& list, const T& val)
{
    const label oldLen = list.size();
    list.resize(oldLen + 1);

    for (label i = oldLen; i > 0; --i)
    {
        list[i] = std::move(list[i-1]);
    }

    list[0] = val;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::boundBox Foam::PDRblock::bounds
(
    const scalarList& x,
    const scalarList& y,
    const scalarList& z
)
{
    return boundBox
    (
        point(x.first(), y.first(), z.first()),
        point(x.last(),  y.last(),  z.last())
    );
}


Foam::Vector<Foam::gradingDescriptors>
Foam::PDRblock::grading(const Vector<gridControl>& ctrl)
{
    return Vector<gradingDescriptors>
    (
        ctrl.x().grading(),
        ctrl.y().grading(),
        ctrl.z().grading()
    );
}

Foam::labelVector
Foam::PDRblock::sizes(const Vector<gridControl>& ctrl)
{
    return labelVector
    (
        ctrl.x().nCells(),
        ctrl.y().nCells(),
        ctrl.z().nCells()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::PDRblock::gridControl::nCells() const
{
    label nTotal = 0;
    for (const label nDiv : divisions_)
    {
        nTotal += nDiv;
    }

    return nTotal;
}


Foam::gradingDescriptors Foam::PDRblock::gridControl::grading() const
{
    // Begin/end nodes for each segment
    const scalarList& knots = *this;

    gradingDescriptors gds(divisions_.size());

    forAll(gds, i)
    {
        //- Construct from components
        gds[i] = gradingDescriptor
        (
            knots[i+1] - knots[i],  // blockFraction from delta
            divisions_[i],          // nDivFraction  from nDivs
            expansion_[i]
        );
    }

    gds.normalise();

    return gds;
}


void Foam::PDRblock::gridControl::resize(label len)
{
    // Begin/end nodes for each segment
    scalarList& knots = *this;

    knots.resize(len, Zero);

    len = Foam::max(0, len-1);

    divisions_.resize(len, Zero);
    expansion_.resize(len, Zero);
}


void Foam::PDRblock::gridControl::append
(
    const scalar p,
    const label nDiv,
    scalar expRatio
)
{
    // Begin/end nodes for each segment
    scalarList& knots = *this;

    // Is monotonic?
    if (knots.size() && (p <= knots.last()))
    {
        WarningInFunction
            << "Cannot append point " << p
            << " which is <= last value " << knots.last() << endl;
        return;
    }

    if (nDiv < 1)
    {
        WarningInFunction
            << "Negative or zero divisions " << nDiv << endl;
        return;
    }

    // Correct expansion ratios - negative is the same as inverse.
    if (expRatio < 0)
    {
        expRatio = 1.0/(-expRatio);
    }
    else if (equal(expRatio, 0))
    {
        expRatio = 1;
    }

    // Now append (push_back)
    knots.append(p);
    divisions_.append(nDiv);
    expansion_.append(expRatio);
}


void Foam::PDRblock::gridControl::prepend
(
    const scalar p,
    const label nDiv,
    scalar expRatio
)
{
    // Begin/end nodes for each segment
    scalarList& knots = static_cast<scalarList&>(*this);

    // Is monotonic?
    if (knots.size() && (p >= knots.first()))
    {
        WarningInFunction
            << "Cannot prepend point " << p
            << " which is >= first value " << knots.first() << endl;
        return;
    }

    if (nDiv < 1)
    {
        WarningInFunction
            << "Negative or zero divisions " << nDiv << endl;
        return;
    }

    // Correct expansion ratios - negative is the same as inverse.
    if (expRatio < 0)
    {
        expRatio = 1.0/(-expRatio);
    }
    else if (equal(expRatio, 0))
    {
        expRatio = 1;
    }

    // Now prepend (push_front)
    prependList(knots, p);
    prependList(divisions_, nDiv);
    prependList(expansion_, expRatio);
}


void Foam::PDRblock::gridControl::writeDict
(
    Ostream& os,
    const direction cmpt
) const
{
    if (cmpt < vector::nComponents)
    {
        os.beginBlock(vector::componentNames[cmpt]);
    }


    const scalarList& knots = *this;

    os  << indent << "points  " << flatOutput(knots);
    os.endEntry();

    os  << indent << "nCells  " << flatOutput(divisions_);
    os.endEntry();

    os  << indent << "ratios  " << flatOutput(expansion_);
    os.endEntry();

    if (cmpt < vector::nComponents)
    {
        os.endBlock();
    }
    os << nl;
}


void Foam::PDRblock::location::reset
(
    const scalar low,
    const scalar upp,
    const label nCells
)
{
    scalarList& grid = *this;

    grid.resize_nocopy(nCells+1);

    grid.front() = low;
    grid.back() = upp;

    const scalar span = (upp - low);

    for (label pointi = 1; pointi < nCells; ++pointi)
    {
        grid[pointi] = low + (pointi*span)/nCells;
    }
}


Foam::scalarMinMax Foam::PDRblock::location::edgeLimits() const
{
    scalarMinMax limits;

    for (label edgei = 0; edgei < this->nCells(); ++edgei)
    {
        limits.add(width(edgei));
    }

    return limits;
}


Foam::label Foam::PDRblock::location::findCell(const scalar p) const
{
    if (scalarList::empty() || p < first() || p > last())
    {
        return -1;
    }
    else if (equal(p, first()))
    {
        return 0;
    }
    else if (equal(p, last()))
    {
        return nCells()-1;
    }
    else if (p < first() || p > last())
    {
        // The point is out-of-bounds
        return -1;
    }

    // Binary search, finds lower index and thus corresponds to the
    // cell in which the point is found
    return findLower(*this, p);
}


Foam::label Foam::PDRblock::location::findIndex
(
    const scalar p,
    const scalar tol
) const
{
    if (scalarList::empty())
    {
        return -1;
    }
    else if (Foam::mag(p - first()) <= tol)
    {
        return 0;
    }
    else if (Foam::mag(p - last()) <= tol)
    {
        return scalarList::size()-1;
    }
    else if (p < first() || p > last())
    {
        // The point is out-of-bounds
        return -1;
    }

    // Linear search
    label i = 0;
    scalar delta = GREAT;

    for (const scalar& val : *this)
    {
        const scalar diff = mag(p - val);

        if (diff <= tol)
        {
            return i;
        }
        else if (delta < diff)
        {
            // Moving further away
            break;
        }

        delta = diff;
        ++i;
    }

    // This point is within bounds, but not near a grid-point
    return -2;
}


// ************************************************************************* //
