/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "FixedList.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T, unsigned N>
std::streamsize Foam::FixedList<T, N>::byteSize()
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return FixedList<T, N>::size_bytes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::find(const T& val) const
{
    const auto iter = std::find(this->cbegin(), this->cend(), val);
    return (iter != this->cend() ? label(iter - this->cbegin()) : label(-1));
}


template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::find
(
    const T& val,
    label pos,
    label len
) const
{
    if (pos >= 0 && pos < label(N))
    {
        // Change sub-length to (one-past) end position
        // len == -1 (like std::string::npos) - search until end

        if (len > 0) len += pos;
        if (len < 0 || len > label(N))
        {
            len = label(N);
        }

        const auto iter = std::find
        (
            (this->cbegin() + pos),
            (this->cbegin() + len),
            val
        );

        if (iter != (this->cbegin() + len))
        {
            return label(iter - this->cbegin());
        }
    }

    return -1;
}


template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::rfind(const T& val, label pos) const
{
    // pos == -1 (like std::string::npos) - search from end

    if (pos < 0 || pos >= label(N))
    {
        pos = label(N)-1;
    }

    while (pos >= 0)
    {
        if (this->v_[pos] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::moveFirst(const label i)
{
    checkIndex(i);

    for (label lower = 0; lower < i; ++lower)
    {
        Foam::Swap(v_[lower], v_[i]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::moveLast(const label i)
{
    checkIndex(i);

    for (label upper = label(N-1); upper > i; --upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::swapFirst(const label i)
{
    checkIndex(i);

    if (i > 0)
    {
        Foam::Swap(v_[0], v_[i]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::swapLast(const label i)
{
    checkIndex(i);

    const label upper = label(N-1);

    if (i < upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator==(const FixedList<T, N>& list) const
{
    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    return
    (
        // List sizes are identical by definition (template parameter)
        std::equal(this->cbegin(), this->cend(), list.cbegin())
    );
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator<(const FixedList<T, N>& list) const
{
    // List sizes are identical by definition (template parameter)

    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    return std::lexicographical_compare
    (
        this->cbegin(), this->cend(),
        list.cbegin(), list.cend()
    );
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator!=(const FixedList<T, N>& list) const
{
    return !operator==(list);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator>(const FixedList<T, N>& list) const
{
    return list.operator<(*this);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator<=(const FixedList<T, N>& list) const
{
    return !list.operator<(*this);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator>=(const FixedList<T, N>& list) const
{
    return !operator<(list);
}


// ************************************************************************* //
