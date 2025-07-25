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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Location

inline bool Foam::PDRblock::location::good() const noexcept
{
    return (scalarList::size() > 1);
}


inline Foam::label Foam::PDRblock::location::nCells() const noexcept
{
    return (scalarList::size()-1);
}


inline Foam::label Foam::PDRblock::location::nPoints() const noexcept
{
    return scalarList::size();
}


inline bool Foam::PDRblock::location::contains(const scalar p) const
{
    return (scalarList::size() > 1 && front() <= p && p <= back());
}


inline const Foam::scalar& Foam::PDRblock::location::min() const
{
    return scalarList::empty() ? pTraits<scalar>::rootMax : front();
}


inline const Foam::scalar& Foam::PDRblock::location::max() const
{
    return scalarList::empty() ? pTraits<scalar>::rootMin : back();
}


inline Foam::scalar Foam::PDRblock::location::centre() const
{
    return scalarList::empty() ? 0 : (0.5*front() + 0.5*back());
}


inline Foam::scalar Foam::PDRblock::location::length() const
{
    return scalarList::empty() ? 0 : mag(back() - front());
}


inline void Foam::PDRblock::location::checkIndex(const label i) const
{
    if (i < 0 || i >= nCells())
    {
        FatalErrorInFunction
            << "The index " << i
            << " is out of range [0," << nCells() << ']' << nl
            << abort(FatalError);
    }
}


inline Foam::scalar Foam::PDRblock::location::width(const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return (operator[](i+1) - operator[](i));
}


inline Foam::scalar Foam::PDRblock::location::C(const label i) const
{
    if (i == -1)
    {
        #ifdef FULLDEBUG
        checkIndex(0);
        #endif

        // "Halo" centre [-1] == x0 - 1/2 (x1 - x0)
        return front() - 0.5*(operator[](1) - front());
    }
    else if (i > 1 && i == scalarList::size()-1)
    {
        // "Halo" centre [nCells] == xN + 1/2 (xN - xN1)
        return back() + 0.5*(back() - operator[](scalarList::size()-2));
    }

    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return 0.5*(operator[](i+1) + operator[](i));
}


inline const Foam::scalar&
Foam::PDRblock::location::clip(const scalar& val) const
{
    if (scalarList::size())
    {
        if (val < front())
        {
            return front();
        }
        else if (back() < val)
        {
            return back();
        }
    }

    return val; // Pass-through
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::PDRblock::dx(const label i) const
{
    return grid_.x().width(i);
}


inline Foam::scalar Foam::PDRblock::dx(const labelVector& ijk) const
{
    return grid_.x().width(ijk.x());
}


inline Foam::scalar Foam::PDRblock::dy(const label j) const
{
    return grid_.y().width(j);
}


inline Foam::scalar Foam::PDRblock::dy(const labelVector& ijk) const
{
    return grid_.y().width(ijk.y());
}


inline Foam::scalar Foam::PDRblock::dz(const label k) const
{
    return grid_.z().width(k);
}


inline Foam::scalar Foam::PDRblock::dz(const labelVector& ijk) const
{
    return grid_.z().width(ijk.z());
}


inline Foam::vector Foam::PDRblock::span
(
    const label i,
    const label j,
    const label k
) const
{
    return vector(dx(i), dy(j), dz(k));
}


inline Foam::vector Foam::PDRblock::span(const labelVector& ijk) const
{
    return vector(dx(ijk), dy(ijk), dz(ijk));
}


inline Foam::point Foam::PDRblock::grid
(
    const label i,
    const label j,
    const label k
) const
{
    return point(grid_.x()[i], grid_.y()[j], grid_.z()[k]);
}


inline Foam::point Foam::PDRblock::grid(const labelVector& ijk) const
{
    return
        point
        (
            grid_.x()[ijk.x()],
            grid_.y()[ijk.y()],
            grid_.z()[ijk.z()]
        );
}


inline Foam::point Foam::PDRblock::C
(
    const label i,
    const label j,
    const label k
) const
{
    return point(grid_.x().C(i), grid_.y().C(j), grid_.z().C(k));
}


inline Foam::point Foam::PDRblock::C(const labelVector& ijk) const
{
    return
        point
        (
            grid_.x().C(ijk.x()),
            grid_.y().C(ijk.y()),
            grid_.z().C(ijk.z())
        );
}


inline Foam::scalar Foam::PDRblock::V
(
    const label i,
    const label j,
    const label k
) const
{
    return dx(i)*dy(j)*dz(k);
}


inline Foam::scalar Foam::PDRblock::V(const labelVector& ijk) const
{
    return dx(ijk.x())*dy(ijk.y())*dz(ijk.z());
}


inline Foam::scalar Foam::PDRblock::width
(
    const label i,
    const label j,
    const label k
) const
{
    return Foam::cbrt(V(i, j, k));
}


inline Foam::scalar Foam::PDRblock::width(const labelVector& ijk) const
{
    return Foam::cbrt(V(ijk));
}


// ************************************************************************* //
