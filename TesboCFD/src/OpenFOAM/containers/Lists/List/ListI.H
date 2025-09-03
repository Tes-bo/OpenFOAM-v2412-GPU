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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::doAlloc()
{
    if (this->size_ > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        this->v_ = new T[this->size_];
    }
}


template<class T>
inline void Foam::List<T>::reAlloc(const label len)
{
    if (this->size_ != len)
    {
        clear();
        this->size_ = len;
        doAlloc();
    }
}


template<class T>
template<class ListType>
inline void Foam::List<T>::copyList(const ListType& list)
{
    // NB: operator[] for list read access (eg, an indirect list)
    // cannot necessarily replace with std::copy

    const label len = this->size_;

    auto iter = this->v_;

    for (label i = 0; i < len; (void)++i, (void)++iter)
    {
        *iter = list[i];
    }
}


template<class T>
template<class ListType, class ListIndices>
inline void Foam::List<T>::copyList
(
    const ListType& list,
    const ListIndices& indices
)
{
    // NB: operator[] for list read access (eg, an indirect list)
    // cannot necessarily replace with std::copy

    const label len = this->size_;

    auto iter = this->v_;

    for (label i = 0; i < len; (void)++i, (void)++iter)
    {
        *iter = list[indices[i]];
    }
}


template<class T>
template<class InputIterator>
inline Foam::List<T>::List
(
    InputIterator firstIter,
    InputIterator lastIter,  // (unused)
    const label len
)
:
    UList<T>(nullptr, len)
{
    if (this->size_)
    {
        doAlloc();

        // Like std::copy() or std::copy_n()
        // but without any requirements on the iterator category

        for (label i = 0; i < len; (void)++i, (void)++firstIter)
        {
            this->v_[i] = *firstIter;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::List<T>::List() noexcept
{}


template<class T>
inline Foam::autoPtr<Foam::List<T>> Foam::List<T>::clone() const
{
    return autoPtr<List<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::clear()
{
    if (this->v_)
    {
        delete[] this->v_;
        this->v_ = nullptr;
    }
    this->size_ = 0;
}


namespace Foam
{
    // Template specialization for bool. Fills new entries with false
    template<>
    inline void List<bool>::resize(const label newLen)
    {
        this->resize(newLen, false);
    }
}


template<class T>
inline void Foam::List<T>::resize(const label len)
{
    this->doResize(len);
}


template<class T>
inline void Foam::List<T>::resize_fill(const label len, const T& val)
{
    this->reAlloc(len);
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::resize_nocopy(const label len)
{
    this->reAlloc(len);
}


template<class T>
inline void Foam::List<T>::resize_unsafe(const label len) noexcept
{
    UList<T>::setAddressableSize(len);
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    label n = this->size();

    if (i >= n)
    {
        if (!n) n = 1;  // Avoid dead-lock when starting from zero-sized

        do
        {
            n *= 2;
        }
        while (i >= n);

        resize(n);
    }

    return UList<T>::operator[](i);
}


template<class T>
template<class... Args>
inline T& Foam::List<T>::emplace_back(Args&&... args)
{
    // This could/should be better with inplace construction
    // (as per std::vector), but currently lacking the methods for that
    // so resize and move assign

    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = T(std::forward<Args>(args)...);
    return UList<T>::operator[](idx);
}


template<class T>
inline void Foam::List<T>::push_back(const T& val)
{
    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = val;  // copy element
}


template<class T>
inline void Foam::List<T>::push_back(T&& val)
{
    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = std::move(val);  // move assign element
}


template<class T>
inline void Foam::List<T>::push_back(const UList<T>& list)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self" << abort(FatalError);
    }

    const label idx = this->size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T>
template<class Addr>
inline void Foam::List<T>::push_back(const IndirectListBase<T, Addr>& list)
{
    // Note: push_back will still work even if the indirect list
    // actually references *this, since its source elements will not
    // overlap the new destinations.

    const label idx = this->size();
    const label n = list.size();

    resize(idx + n);

    auto iter = this->begin(idx);

    for (label i = 0; i < n; (void)++i, (void)++iter)
    {
        *iter = list[i];  // copy element
    }
}


template<class T>
inline Foam::label Foam::List<T>::push_uniq(const T& val)
{
    if (this->contains(val))
    {
        return 0;
    }
    else
    {
        this->push_back(val);
        return 1;  // Increased list length by one
    }
}


template<class T>
inline void Foam::List<T>::pop_back(label n)
{
    if (n >= this->size())
    {
        this->clear();
    }
    else if (n > 0)
    {
        resize(this->size() - n);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::operator=(const Foam::zero)
{
    UList<T>::operator=(Foam::zero{});
}


// ************************************************************************* //
