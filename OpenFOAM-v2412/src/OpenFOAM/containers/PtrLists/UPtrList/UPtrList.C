/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "UPtrList.H"
#include "PtrList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::UPtrList<T>::UPtrList(PtrList<T>& list)
:
    ptrs_(list.ptrs_) // shallow copy (via const reference)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::UPtrList<T>::squeezeNull()
{
    const label len = this->size();
    label newLen = 0;

    for (label i=0; i < len; ++i)
    {
        T* ptr = ptrs_[i];
        if (ptr)
        {
            if (i != newLen)
            {
                ptrs_[newLen] = ptr;
                ptrs_[i] = nullptr;
            }
            ++newLen;
        }
    }

    return newLen;
}


template<class T>
void Foam::UPtrList<T>::trimTrailingNull()
{
    label newLen = this->size();

    for (label i = newLen-1; i >= 0 && !ptrs_[i]; --i)
    {
        --newLen;
    }

    // Or mutable?
    // const_cast<Detail::PtrListDetail<T>&>(ptrs_).setAddressableSize(newLen);

    ptrs_.setAddressableSize(newLen);
}


template<class T>
void Foam::UPtrList<T>::reorder(const labelUList& oldToNew, const bool check)
{
    const label len = this->size();

    if (oldToNew.size() != len)
    {
        FatalErrorInFunction
            << "Size of map (" << oldToNew.size()
            << ") not equal to list size (" << len
            << ") for type " << typeid(T).name() << nl
            << abort(FatalError);
    }

    Detail::PtrListDetail<T> newList(len);

    for (label i=0; i<len; ++i)
    {
        const label newIdx = oldToNew[i];

        if (newIdx < 0 || newIdx >= len)
        {
            FatalErrorInFunction
                << "Illegal index " << newIdx << nl
                << "Valid indices are [0," << len << ") for type "
                << typeid(T).name() << nl
                << abort(FatalError);
        }

        if (newList[newIdx])
        {
            FatalErrorInFunction
                << "reorder map is not unique; element " << newIdx
                << " already used for type " << typeid(T).name()
                << abort(FatalError);
        }
        newList[newIdx] = ptrs_[i];
    }

    // Verify all pointers were indeed set
    if (check)
    {
        newList.checkNonNull();
    }

    ptrs_.transfer(newList);
}


template<class T>
void Foam::UPtrList<T>::sortOrder(const labelUList& order, const bool check)
{
    const label len = this->size();

    if (order.size() != len)
    {
        FatalErrorInFunction
            << "Size of map (" << order.size()
            << ") not equal to list size (" << len
            << ") for type " << typeid(T).name() << nl
            << abort(FatalError);
    }

    Detail::PtrListDetail<T> newList(len);
    Detail::PtrListDetail<T> guard(len);

    for (label i=0; i<len; ++i)
    {
        const label oldIdx = order[i];

        if (oldIdx < 0 || oldIdx >= len)
        {
            FatalErrorInFunction
                << "Illegal index " << oldIdx << nl
                << "Valid indices are [0," << len << ") for type "
                << typeid(T).name() << nl
                << abort(FatalError);
        }

        if (guard[oldIdx])
        {
            FatalErrorInFunction
                << "order map is not unique; element " << oldIdx
                << " already used for type " << typeid(T).name()
                << abort(FatalError);
        }

        guard[oldIdx] = ptrs_[oldIdx];
        newList[i] = ptrs_[oldIdx];
    }

    // Verify that all pointers were indeed set
    if (check)
    {
        newList.checkNonNull();
    }

    ptrs_.transfer(newList);
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const UPtrList<T>& list)
{
    return list.ptrs_.write(os);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::sort(UPtrList<T>& list)
{
    std::stable_sort
    (
        list.begin_ptr(),
        list.end_ptr(),
        // Compare less, with nullptr protect and sort nullptr to end
        [](const T* const a, const T* const b) -> bool
        {
            return (a && b) ? (*a < *b) : !b;
        }
    );
}


template<class T, class Compare>
void Foam::sort(UPtrList<T>& list, const Compare& comp)
{
    std::stable_sort
    (
        list.begin_ptr(),
        list.end_ptr(),
        typename UPtrList<T>::template value_compare<Compare>(comp)
    );
}


// ************************************************************************* //
