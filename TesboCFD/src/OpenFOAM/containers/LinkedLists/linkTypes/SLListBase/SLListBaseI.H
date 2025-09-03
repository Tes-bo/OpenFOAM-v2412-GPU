/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "error.H"
#include "nullObject.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class IteratorType>
inline const IteratorType& Foam::SLListBase::iterator_end()
{
    return *reinterpret_cast<const IteratorType*>(nullObjectPtr);
}


template<class IteratorType>
inline IteratorType Foam::SLListBase::iterator_first() const
{
    SLListBase* list = const_cast<SLListBase*>(this);

    if (size())
    {
        return IteratorType(list, const_cast<SLListBase::link*>(last_->next_));
    }

    // Return an end iterator
    return IteratorType(list, nullptr);
}


template<class IteratorType>
inline IteratorType Foam::SLListBase::iterator_last() const
{
    SLListBase* list = const_cast<SLListBase*>(this);

    if (size())
    {
        return IteratorType(list, const_cast<SLListBase::link*>(last_));
    }

    // Return an end iterator
    return IteratorType(list, nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::SLListBase::link*
Foam::SLListBase::front()
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline const Foam::SLListBase::link*
Foam::SLListBase::front() const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline Foam::SLListBase::link*
Foam::SLListBase::back()
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const Foam::SLListBase::link*
Foam::SLListBase::back() const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void Foam::SLListBase::clear()
{
    last_ = nullptr;
    size_ = 0;
}


inline void Foam::SLListBase::swap(SLListBase& lst)
{
    if (this == &lst)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(last_, lst.last_);
    std::swap(size_, lst.size_);
}


inline void Foam::SLListBase::transfer(SLListBase& lst)
{
    if (this == &lst)
    {
        return;  // Self-assignment is a no-op
    }

    last_ = lst.last_;
    size_ = lst.size_;

    lst.clear();
}


inline Foam::SLListBase::link* Foam::SLListBase::remove
(
    SLListBase::iterator& iter
)
{
    return remove(iter.node_);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline Foam::SLListBase::iterator::iterator
(
    SLListBase* list,
    SLListBase::link* item
)
:
    node_(item),
    list_(list),
    copy_()
{
    if (node_ != nullptr)
    {
        copy_ = *node_;
    }
}


inline Foam::SLListBase::link*
Foam::SLListBase::iterator::get_node() const noexcept
{
    return node_;
}


inline bool Foam::SLListBase::iterator::good() const noexcept
{
    return (node_ != nullptr);
}


inline void Foam::SLListBase::iterator::next()
{
    if (list_)
    {
        if (node_ == list_->last_ || list_->last_ == nullptr)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = copy_.next_;
            copy_ = *node_;
        }
    }
}


inline void Foam::SLListBase::iterator::operator=(const iterator& iter)
{
    node_ = iter.node_;
    list_ = iter.list_;
    copy_ = iter.copy_;
}


inline bool Foam::SLListBase::iterator::operator==(const iterator& iter) const
{
    return node_ == iter.node_;
}


inline bool Foam::SLListBase::iterator::operator!=(const iterator& iter) const
{
    return node_ != iter.node_;
}


inline Foam::SLListBase::iterator
Foam::SLListBase::begin()
{
    if (size())
    {
        return iterator_first<iterator>();
    }

    return end();
}


inline const Foam::SLListBase::iterator&
Foam::SLListBase::end()
{
    return iterator_end<SLListBase::iterator>();
}


inline const Foam::SLListBase::const_iterator&
Foam::SLListBase::cend() const
{
    return iterator_end<SLListBase::const_iterator>();
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline Foam::SLListBase::const_iterator::const_iterator
(
    const SLListBase* list,
    const SLListBase::link* item
)
:
    node_(item),
    list_(list)
{}


inline Foam::SLListBase::const_iterator::const_iterator
(
    const SLListBase::iterator& iter
)
:
    node_(iter.node_),
    list_(iter.list_)
{}


inline const Foam::SLListBase::link*
Foam::SLListBase::const_iterator::get_node() const noexcept
{
    return node_;
}


inline bool Foam::SLListBase::const_iterator::good() const noexcept
{
    return (node_ != nullptr);
}


inline void Foam::SLListBase::const_iterator::next()
{
    if (list_)
    {
        if (node_ == list_->last_)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = node_->next_;
        }
    }
}


inline bool Foam::SLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return node_ == iter.node_;
}


inline bool Foam::SLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return node_ != iter.node_;
}


inline Foam::SLListBase::const_iterator
Foam::SLListBase::cbegin() const
{
    if (size())
    {
        return iterator_first<const_iterator>();
    }

    return cend();
}


// ************************************************************************* //
