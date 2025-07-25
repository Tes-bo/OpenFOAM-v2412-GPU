/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::UPtrList<T>::setAddressableSize(const label n) noexcept
{
    ptrs_.setAddressableSize(n);
}


template<class T>
inline Foam::label Foam::UPtrList<T>::find_next(label pos) const
{
    return ptrs_.find_next(pos);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::UPtrList<T>::UPtrList(const label len)
:
    ptrs_(len)
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(Detail::PtrListDetail<T>&& ptrs) noexcept
:
    ptrs_(std::move(ptrs))
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(const UPtrList<T>& list)
:
    ptrs_(list.ptrs_)
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(UPtrList<T>&& list) noexcept
:
    ptrs_(std::move(list.ptrs_))
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(UPtrList<T>& list, bool reuse)
:
    ptrs_(list.ptrs_, reuse)
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(const UList<T*>& list)
:
    ptrs_(list)
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(UList<T>& list)
:
    ptrs_(list.size())
{
    const label len = ptrs_.size();

    for (label i = 0; i < len; ++i)
    {
        ptrs_[i] = &(list[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::UPtrList<T>::empty() const noexcept
{
    return ptrs_.empty();
}


template<class T>
inline Foam::label Foam::UPtrList<T>::size() const noexcept
{
    return ptrs_.size();
}


template<class T>
inline Foam::label Foam::UPtrList<T>::capacity() const noexcept
{
    return ptrs_.capacity();
}


template<class T>
inline Foam::label Foam::UPtrList<T>::count_nonnull() const noexcept
{
    return ptrs_.count_nonnull();
}


template<class T>
inline const T* Foam::UPtrList<T>::test(const label i) const
{
    return ptrs_.get(i);
}


template<class T>
inline const T* Foam::UPtrList<T>::get(const label i) const
{
    return ptrs_.get(i);
}


template<class T>
inline T* Foam::UPtrList<T>::get(const label i)
{
    return ptrs_.get(i);
}


template<class T>
inline const T& Foam::UPtrList<T>::at(const label i) const
{
    const T* ptr = ptrs_.get(i);

    if (!ptr)
    {
        FatalErrorInFunction
            << "Cannot dereference nullptr at index " << i
            << " in range [0," << size() << ")\n"
            << abort(FatalError);
    }

    return *ptr;
}


template<class T>
inline T& Foam::UPtrList<T>::at(const label i)
{
    T* ptr = ptrs_.get(i);

    if (!ptr)
    {
        FatalErrorInFunction
            << "Cannot dereference nullptr at index " << i
            << " in range [0," << size() << ")\n"
            << abort(FatalError);
    }

    return *ptr;
}


template<class T>
inline T* Foam::UPtrList<T>::set(const label i, T* ptr)
{
    T* old = ptrs_[i];
    if (old == ptr)
    {
        return nullptr;  // Content did not change
    }
    ptrs_[i] = ptr;
    return old;
}


template<class T>
inline void Foam::UPtrList<T>::clear()
{
    ptrs_.clear();
}


template<class T>
inline void Foam::UPtrList<T>::free()
{
    ptrs_ = nullptr;
}


template<class T>
inline void Foam::UPtrList<T>::swap(UPtrList<T>& list) noexcept
{
    ptrs_.swap(list.ptrs_);
}


template<class T>
inline void Foam::UPtrList<T>::transfer(UPtrList<T>& list)
{
    ptrs_.transfer(list.ptrs_);
}


template<class T>
inline T& Foam::UPtrList<T>::front()
{
    return this->at(0);
}


template<class T>
inline const T& Foam::UPtrList<T>::front() const
{
    return this->at(0);
}


template<class T>
inline T& Foam::UPtrList<T>::back()
{
    return this->at(this->size()-1);
}


template<class T>
inline const T& Foam::UPtrList<T>::back() const
{
    return this->at(this->size()-1);
}


template<class T>
inline void Foam::UPtrList<T>::resize(const label newLen)
{
    ptrs_.resize(newLen);
}


template<class T>
inline void Foam::UPtrList<T>::resize_null(const label newLen)
{
    ptrs_.resize_null(newLen);
}


template<class T>
inline void Foam::UPtrList<T>::push_back(T* ptr)
{
    ptrs_.push_back(ptr);
}


template<class T>
inline void Foam::UPtrList<T>::push_back(UPtrList<T>&& other)
{
    ptrs_.push_back(other.ptrs_);
    other.ptrs_.clear();
}


template<class T>
inline void Foam::UPtrList<T>::checkNonNull() const
{
    ptrs_.checkNonNull();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::UPtrList<T>::operator[](const label i) const
{
    return this->at(i);
}


template<class T>
inline T& Foam::UPtrList<T>::operator[](const label i)
{
    return this->at(i);
}


// * * * * * * * * * * * * * * * iterator base * * * * * * * * * * * * * * * //

template<class T>
template<bool Const>
inline constexpr
Foam::UPtrList<T>::Iterator<Const>::Iterator() noexcept
:
    list_(nullptr),
    pos_(-1)
{}


template<class T>
template<bool Const>
inline Foam::UPtrList<T>::Iterator<Const>::Iterator
(
    list_type* list
)
:
    list_(list),
    pos_(-1)
{
    if (list_)
    {
        pos_ = list_->find_next(-1);
    }
}


template<class T>
template<bool Const>
inline void Foam::UPtrList<T>::Iterator<Const>::increment()
{
    //TDB: extra safety?  if (iter.good())
    {
        pos_ = list_->find_next(pos_);
    }
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline T* Foam::UPtrList<T>::iterator::get() const
{
    return this->good() ? this->list_->get(this->pos_) : nullptr;
}


template<class T>
inline T& Foam::UPtrList<T>::iterator::val() const
{
    return this->list_->at(this->pos_);
}


template<class T>
inline typename Foam::UPtrList<T>::iterator&
Foam::UPtrList<T>::iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator++(int)
{
    iterator iter(*this);
    this->increment();
    return iter;
}


// * * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * //

template<class T>
inline const T* Foam::UPtrList<T>::const_iterator::get() const
{
    return this->good() ? this->list_->get(this->pos_) : nullptr;
}


template<class T>
inline const T& Foam::UPtrList<T>::const_iterator::val() const
{
    return this->list_->at(this->pos_);
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator&
Foam::UPtrList<T>::const_iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::const_iterator::operator++(int)
{
    const_iterator iter(*this);
    this->increment();
    return iter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::begin()
{
    return iterator(Iterator<false>(this));
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::begin() const
{
    return const_iterator(Iterator<true>(this));
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::cbegin() const
{
    return const_iterator(Iterator<true>(this));
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::end() noexcept
{
    return iterator();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::end() const noexcept
{
    return UPtrList<T>::const_iterator();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::cend() const noexcept
{
    return UPtrList<T>::const_iterator();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::UPtrList<T>::operator=(const UPtrList<T>& list)
{
    ptrs_ = list.ptrs_;  // shallow copy
}


template<class T>
inline void Foam::UPtrList<T>::operator=(UPtrList<T>&& list)
{
    ptrs_.transfer(list.ptrs_);
}


// ************************************************************************* //
