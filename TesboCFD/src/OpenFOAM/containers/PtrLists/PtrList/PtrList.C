/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "PtrList.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::PtrList(PtrList<T>& list, bool reuse)
:
    UPtrList<T>(list, reuse)
{
    if (!reuse)
    {
        // This works like an inplace clone method
        const label len = this->size();

        for (label i=0; i<len; ++i)
        {
            this->ptrs_[i] = (list[i]).clone().ptr();
        }
    }
}


template<class T>
Foam::PtrList<T>::PtrList(const SLPtrList<T>& list)
:
    UPtrList<T>(list.size())
{
    if (list.size())
    {
        label i = 0;
        for (auto iter = list.cbegin(); iter != list.cend(); ++iter)
        {
            this->ptrs_[i++] = (*iter).clone().ptr();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::~PtrList()
{
    (this->ptrs_).free();  // Free old pointers
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
template<class... Args>
Foam::PtrList<T> Foam::PtrList<T>::clone(Args&&... args) const
{
    const label len = this->size();

    PtrList<T> cloned(len);

    for (label i=0; i<len; ++i)
    {
        const T* ptr = this->ptrs_[i];

        if (ptr)
        {
            cloned.ptrs_[i] = ptr->clone(std::forward<Args>(args)...).ptr();
        }
    }

    return cloned;
}


template<class T>
void Foam::PtrList<T>::resize(const label newLen)
{
    const label oldLen = this->size();

    if (newLen <= 0)
    {
        clear();
    }
    else if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i = newLen; i < oldLen; ++i)
        {
            delete this->ptrs_[i];
            this->ptrs_[i] = nullptr;
        }

        // Any new elements are initialized to nullptr.
        (this->ptrs_).resize(newLen);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::operator=(const PtrList<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    const label oldLen = this->size();
    const label newLen = list.size();

    // Truncate (frees old pointers) or extend the length
    resize(newLen);

    if (newLen < oldLen)
    {
        // Copy values for existing entries
        for (label i=0; i<newLen; ++i)
        {
            (*this)[i] = list[i];
        }
    }
    else
    {
        // Copy values for existing entries
        for (label i=0; i<oldLen; ++i)
        {
            (*this)[i] = list[i];
        }

        // Clone pointers for new entries
        for (label i=oldLen; i<newLen; ++i)
        {
            this->ptrs_[i] = (list[i]).clone().ptr();
        }
    }
}


// ************************************************************************* //
