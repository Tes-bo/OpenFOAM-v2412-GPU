/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "ListOps.H"
#include "SubList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::CompactListList<T>::enforceSizeSanity()
{
    if (offsets_.size() == 1)
    {
        offsets_.clear();
    }
    if (offsets_.empty())
    {
        values_.clear();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const CompactListList<T>& list
)
:
    offsets_(list.offsets_),
    values_(list.values_)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>&& list
)
:
    offsets_(std::move(list.offsets_)),
    values_(std::move(list.values_))
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>& list,
    bool reuse
)
:
    offsets_(list.offsets_, reuse),
    values_(list.values_, reuse)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals
)
:
    offsets_(mRows+1, Foam::zero{}),
    values_(nVals)
{
    // Optionally: enforceSizeSanity();
}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals,
    const Foam::zero
)
:
    offsets_(mRows+1, Foam::zero{}),
    values_(nVals, Foam::zero{})
{
    // Optionally: enforceSizeSanity();
}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals,
    const T& val
)
:
    offsets_(mRows+1, Foam::zero{}),
    values_(nVals, val)
{
    // Optionally: enforceSizeSanity();
}


template<class T>
inline Foam::autoPtr<Foam::CompactListList<T>>
Foam::CompactListList<T>::clone() const
{
    return autoPtr<CompactListList<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T* Foam::CompactListList<T>::cdata() const noexcept
{
    return values_.cdata();
}


template<class T>
inline T* Foam::CompactListList<T>::data() noexcept
{
    return values_.data();
}


template<class T>
inline const char* Foam::CompactListList<T>::cdata_bytes() const noexcept
{
    return values_.cdata_bytes();
}


template<class T>
inline char* Foam::CompactListList<T>::data_bytes() noexcept
{
    return values_.data_bytes();
}


template<class T>
inline std::streamsize Foam::CompactListList<T>::size_bytes() const noexcept
{
    return values_.size_bytes();
}


template<class T>
inline bool Foam::CompactListList<T>::empty() const noexcept
{
    // Note: could (should?) also check total size??
    // const label len = (offsets_.size() - 1);
    // return (len < 1) || (*(offsets_.cdata() + len) == 0);
    return (offsets_.size() <= 1);
}


template<class T>
inline bool Foam::CompactListList<T>::single() const noexcept
{
    return (offsets_.size() == 2);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::length() const noexcept
{
    const label len = (offsets_.size() - 1);
    return (len < 1) ? static_cast<label>(0) : len;
}


template<class T>
inline Foam::label Foam::CompactListList<T>::size() const noexcept
{
    const label len = (offsets_.size() - 1);
    return (len < 1) ? static_cast<label>(0) : len;
}


template<class T>
inline Foam::labelList Foam::CompactListList<T>::sizes() const
{
    return localSizes();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::totalSize() const noexcept
{
    return offsets_.empty() ? 0 : *(offsets_.cdata() + offsets_.size()-1);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::maxSize() const
{
    return this->maxNonLocalSize(-1);
}


template<class T>
inline const Foam::labelUList
Foam::CompactListList<T>::localStarts() const
{
    const label len = (offsets_.size() - 1);

    if (len < 1) return labelUList::null();

    return labelList::subList(offsets_, len);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localStart(const label i) const
{
    return offsets_[i];
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localEnd(const label i) const
{
    return offsets_[i+1];
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localSize(const label i) const
{
    return offsets_[i+1] - offsets_[i];
}


template<class T>
inline const Foam::SubList<T>
Foam::CompactListList<T>::localList(const label i) const
{
    return SubList<T>(values_, (offsets_[i+1] - offsets_[i]), offsets_[i]);
}


template<class T>
inline Foam::SubList<T>
Foam::CompactListList<T>::localList(const label i)
{
    return SubList<T>(values_, (offsets_[i+1] - offsets_[i]), offsets_[i]);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::toGlobal
(
    const label rowi,
    const label i
) const
{
    return i + offsets_[rowi];
}


template<class T>
inline Foam::label Foam::CompactListList<T>::toLocal
(
    const label rowi,
    const label i
) const
{
    const label locali = i - offsets_[rowi];

    if (locali < 0 || i >= offsets_[rowi+1])
    {
        FatalErrorInFunction
            << "Index " << i << " does not belong on row "
            << rowi << nl << "Offsets:" << offsets_
            << abort(FatalError);
    }

    return locali;
}


template<class T>
inline Foam::label Foam::CompactListList<T>::findRow(const label i) const
{
    return (i < 0 || i >= totalSize()) ? -1 : findLower(offsets_, i+1);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::whichRow(const label i) const
{
    const label rowi = findRow(i);

    if (rowi < 0)
    {
        FatalErrorInFunction
            << "Index " << i << " outside of range" << nl
            << "Offsets:" << offsets_
            << abort(FatalError);
    }

    return rowi;
}


template<class T>
inline void Foam::CompactListList<T>::clear()
{
    offsets_.clear();
    values_.clear();
}


template<class T>
inline void Foam::CompactListList<T>::resize(const label mRows)
{
    if (mRows == 0)
    {
        // Clear
        offsets_.clear();
        values_.clear();
    }
    else if (mRows < size())
    {
        // Shrink
        offsets_.resize(mRows+1);
        values_.resize(offsets_[mRows]);
    }
    else if (mRows > size())
    {
        // Extend number of rows, each with local size of 0
        const label endOffset = offsets_.empty() ? 0 : offsets_.back();

        offsets_.resize(mRows+1, endOffset);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nVals
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        offsets_.clear();
        values_.clear();
    }
    else
    {
        offsets_.resize(mRows+1, Foam::zero{});
        values_.resize(nVals);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize_nocopy
(
    const label mRows,
    const label nVals
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        offsets_.clear();
        values_.clear();
    }
    else
    {
        offsets_.resize(mRows+1, Foam::zero{});
        values_.resize_nocopy(nVals);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nVals,
    const T& val
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        offsets_.clear();
        values_.clear();
    }
    else
    {
        offsets_.resize(mRows+1, Foam::zero{});
        values_.resize(nVals, val);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::CompactListList<T>::operator=
(
    const CompactListList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    offsets_ = list.offsets_,
    values_ = list.values_;
}


template<class T>
inline void Foam::CompactListList<T>::operator=
(
    CompactListList<T>&& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    offsets_.transfer(list.offsets_);
    values_.transfer(list.values_);
}


template<class T>
inline void Foam::CompactListList<T>::operator=(const T& val)
{
    values_ = val;
}


template<class T>
inline void Foam::CompactListList<T>::operator=(const Foam::zero)
{
    values_ = Foam::zero{};
}


template<class T>
inline const Foam::SubList<T>
Foam::CompactListList<T>::operator[](const label i) const
{
    // return SubList<T>(values_, (offsets_[i+1] - offsets_[i]), offsets_[i]);
    return this->localList(i);
}


template<class T>
inline Foam::SubList<T>
Foam::CompactListList<T>::operator[](const label i)
{
    // return SubList<T>(values_, (offsets_[i+1] - offsets_[i]), offsets_[i]);
    return this->localList(i);
}


template<class T>
inline T& Foam::CompactListList<T>::operator()
(
    const label i,
    const label j
)
{
    return values_[toGlobal(i, j)];
}


template<class T>
inline const T& Foam::CompactListList<T>::operator()
(
    const label i,
    const label j
) const
{
    return values_[toGlobal(i, j)];
}


// ************************************************************************* //
