/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

inline Foam::ijkAddressing::ijkAddressing()
:
    sizes_(0, 0, 0)
{}


inline Foam::ijkAddressing::ijkAddressing(const labelVector& ijk)
:
    sizes_(ijk)
{
    #ifdef FULLDEBUG
    checkSizes();
    #endif
}


inline Foam::ijkAddressing::ijkAddressing
(
    const label ni,
    const label nj,
    const label nk
)
:
    sizes_(ni, nj, nk)
{
    #ifdef FULLDEBUG
    checkSizes();
    #endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::ijkAddressing::empty() const
{
    return (!sizes_.x() || !sizes_.y() || !sizes_.z());
}


inline const Foam::labelVector& Foam::ijkAddressing::sizes() const
{
    return sizes_;
}


inline Foam::labelVector& Foam::ijkAddressing::sizes()
{
    return sizes_;
}


inline Foam::label Foam::ijkAddressing::size() const
{
    // Could also use cmptProduct(sizes_);
    return (sizes_.x() * sizes_.y() * sizes_.z());
}


inline const Foam::label& Foam::ijkAddressing::size
(
    const vector::components cmpt
) const
{
    return sizes_[cmpt];
}


inline void Foam::ijkAddressing::clear()
{
    sizes_ = Zero;
}


inline void Foam::ijkAddressing::reset
(
    const label ni,
    const label nj,
    const label nk
)
{
    sizes_.x() = ni;
    sizes_.y() = nj;
    sizes_.z() = nk;

    #ifdef FULLDEBUG
    checkSizes();
    #endif
}


inline void Foam::ijkAddressing::reset(const labelVector& newSizes)
{
    sizes_ = newSizes;

    #ifdef FULLDEBUG
    checkSizes();
    #endif
}


inline Foam::label Foam::ijkAddressing::index
(
    const label i,
    const label j,
    const label k
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j, k);
    #endif

    return (i + (sizes_.x() * (j + (sizes_.y() * k))));
}


inline Foam::label Foam::ijkAddressing::index(const labelVector& ijk) const
{
    #ifdef FULLDEBUG
    checkIndex(ijk);
    #endif

    return (ijk.x() + (sizes_.x() * (ijk.y() + (sizes_.y() * ijk.z()))));
}


inline Foam::labelVector Foam::ijkAddressing::index(const label idx) const
{
    // No checkIndex

    return labelVector
    (
        idx % (sizes_.x()),
        idx / (sizes_.x()),
        idx / (sizes_.x() * sizes_.y())
    );
}


inline void Foam::ijkAddressing::checkIndex
(
    const label i,
    const label j,
    const label k,
    const bool allowExtra
) const
{
    const label extra = (allowExtra ? 1 : 0);

    if (i < 0 || i >= (sizes_.x() + extra))
    {
        FatalErrorInFunction
            << "The i-index " << i
            << " is out of range [0," << (sizes_.x() + extra) << ']' << nl
            << abort(FatalError);
    }
    if (j < 0 || j >= (sizes_.y() + extra))
    {
        FatalErrorInFunction
            << "The j-index " << j
            << " is out of range [0," << (sizes_.y() + extra) << ']' << nl
            << abort(FatalError);
    }
    if (k < 0 || k >= (sizes_.z() + extra))
    {
        FatalErrorInFunction
            << "The k-index " << k
            << " is out of range [0," << (sizes_.z() + extra) << ']' << nl
            << abort(FatalError);
    }
}


inline void Foam::ijkAddressing::checkIndex
(
    const labelVector& ijk,
    const bool allowExtra
) const
{
    checkIndex(ijk.x(), ijk.y(), ijk.z(), allowExtra);
}


inline void Foam::ijkAddressing::checkSizes() const
{
    if (sizes_.x() < 0)
    {
        FatalErrorInFunction
            << "The i-size is negative" << nl
            << abort(FatalError);
    }
    if (sizes_.y() < 0)
    {
        FatalErrorInFunction
            << "The j-size is negative" << nl
            << abort(FatalError);
    }
    if (sizes_.z() < 0)
    {
        FatalErrorInFunction
            << "The k-size is negative" << nl
            << abort(FatalError);
    }
}


inline void Foam::ijkAddressing::checkSizes(const labelVector& other) const
{
    if (sizes_ != other)
    {
        FatalErrorInFunction
            << "The i-j-k sizes are different. "
            << sizes_ << " vs. " << other << nl
            << abort(FatalError);
    }
}


inline void Foam::ijkAddressing::checkSizes(const label nTotal) const
{
    if (size() != nTotal)
    {
        FatalErrorInFunction
            << "The total size is different. "
            << size() << " vs. " << nTotal << nl
            << abort(FatalError);
    }
}


// ************************************************************************* //
