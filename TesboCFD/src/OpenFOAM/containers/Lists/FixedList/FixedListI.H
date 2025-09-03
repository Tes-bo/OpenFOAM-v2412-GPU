/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "UList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(const T& val)
{
    this->fill(val);
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(const Foam::zero)
{
    this->fill(Foam::zero{});
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(const FixedList<T, N>& list)
{
    // std::execution::unsequenced_policy
    std::copy(list.begin(), list.end(), v_);
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(FixedList<T, N>&& list)
{
    // std::execution::unsequenced_policy
    std::move(list.begin(), list.end(), v_);
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(std::initializer_list<T> list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(const UList<T>& list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
template<unsigned AnyNum>
inline Foam::FixedList<T, N>::FixedList
(
    const FixedList<T, AnyNum>& list,
    const FixedList<label, N>& indices
)
{
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = list[indices[i]];
    }
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList
(
    const UList<T>& list,
    const FixedList<label, N>& indices
)
{
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = list[indices[i]];
    }
}


template<class T, unsigned N>
inline Foam::autoPtr<Foam::FixedList<T, N>>
Foam::FixedList<T, N>::clone() const
{
    return autoPtr<FixedList<T, N>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned N>
inline const T*
Foam::FixedList<T, N>::cdata() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline T*
Foam::FixedList<T, N>::data() noexcept
{
    return v_;
}


template<class T, unsigned N>
inline const char*
Foam::FixedList<T, N>::cdata_bytes() const noexcept
{
    return reinterpret_cast<const char*>(v_);
}


template<class T, unsigned N>
inline char*
Foam::FixedList<T, N>::data_bytes() noexcept
{
    return reinterpret_cast<char*>(v_);
}


template<class T, unsigned N>
inline std::streamsize Foam::FixedList<T, N>::size_bytes() noexcept
{
    return N*sizeof(T);
}


template<class T, unsigned N>
template<unsigned Index>
inline T& Foam::FixedList<T, N>::get() noexcept
{
    static_assert(Index < N, "Address outside FixedList range");
    return v_[Index];
}


template<class T, unsigned N>
template<unsigned Index>
inline const T& Foam::FixedList<T, N>::get() const noexcept
{
    static_assert(Index < N, "Address outside FixedList range");
    return v_[Index];
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::front() noexcept
{
    return v_[0];
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::front() const noexcept
{
    return v_[0];
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::back() noexcept
{
    return v_[N-1];
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::back() const noexcept
{
    return v_[N-1];
}


template<class T, unsigned N>
inline Foam::label Foam::FixedList<T, N>::fcIndex(const label i) const noexcept
{
    return (i == N-1 ? 0 : i+1);
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::fcValue(const label i) const
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::fcValue(const label i)
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned N>
inline Foam::label Foam::FixedList<T, N>::rcIndex(const label i) const noexcept
{
    return (i ? i-1 : N-1);
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::rcValue(const label i) const
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::rcValue(const label i)
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkStart(const label start) const
{
    if (start < 0 || (start && unsigned(start) >= N))
    {
        // Note: always accept start=0, even for zero-sized lists
        FatalErrorInFunction
            << "start " << start << " out of range [0," << N << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkSize(const label size) const
{
    if (unsigned(size) != N)
    {
        FatalErrorInFunction
            << "size " << size << " != " << N
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkIndex(const label i) const
{
    if (i < 0 || unsigned(i) >= N)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0," << N << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::uniform() const
{
    if (empty()) return false;   // <- Compile-time disabled anyhow

    // std::all_of()
    for (unsigned i=1; i<N; ++i)
    {
        if (v_[0] != v_[i])
        {
            return false;
        }
    }

    return true;
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::contains(const T& val) const
{
    const auto iter = std::find(this->cbegin(), this->cend(), val);
    return (iter != this->cend());
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::contains
(
    const T& val,
    label pos,
    label len
) const
{
    return (this->find(val, pos, len) >= 0);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize_fill(const label n, const T& val)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
    this->fill(val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize_nocopy(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::fill(const T& val)
{
    // Usually small enough that parallel execution is pointless...
    std::fill_n(v_, N, val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::fill(const Foam::zero)
{
    // Usually small enough that parallel execution is pointless...
    // Cannot use std::fill (ambiguous conversions for bool, char, etc)

    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = Zero;
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::swap(FixedList<T, N>& other)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    // Essentially std::swap_ranges
    for (unsigned i=0; i<N; ++i)
    {
        Foam::Swap(v_[i], other.v_[i]);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::transfer(FixedList<T, N>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // std::execution::unsequenced_policy
    std::move(list.begin(), list.end(), v_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::operator[](const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const UList<T>& list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(std::initializer_list<T> list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const T& val)
{
    this->fill(val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const Foam::zero)
{
    this->fill(Foam::zero{});
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const FixedList<T, N>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // std::execution::unsequenced_policy
    std::copy(list.begin(), list.end(), v_);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(FixedList<T, N>&& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // std::execution::unsequenced_policy
    std::move(list.begin(), list.end(), v_);
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::iterator
Foam::FixedList<T, N>::begin() noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::begin() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::cbegin() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::iterator
Foam::FixedList<T, N>::end() noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::end() const noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::cend() const noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::reverse_iterator
Foam::FixedList<T, N>::rbegin()
{
    return reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::rbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::crbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::reverse_iterator
Foam::FixedList<T, N>::rend()
{
    return reverse_iterator(begin());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::rend() const
{
    return const_reverse_iterator(begin());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::crend() const
{
    return const_reverse_iterator(begin());
}


// ************************************************************************* //
