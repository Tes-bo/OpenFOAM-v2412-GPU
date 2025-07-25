/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "DLListBase.H"
#include "error.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DLListBase::push_front(DLListBase::link* item)
{
    if (!item)
    {
        return;
    }

    ++size_;

    if (!first_)
    {
        item->prev_ = item;
        item->next_ = item;
        first_ = last_ = item;
    }
    else
    {
        item->prev_ = item;
        item->next_ = first_;
        first_->prev_ = item;
        first_ = item;
    }
}


void Foam::DLListBase::push_back(DLListBase::link* item)
{
    if (!item)
    {
        return;
    }

    ++size_;

    if (!first_)
    {
        item->prev_ = item;
        item->next_ = item;
        first_ = last_ = item;
    }
    else
    {
        last_->next_ = item;
        item->prev_ = last_;
        item->next_ = item;
        last_ = item;
    }
}


bool Foam::DLListBase::swapUp(DLListBase::link* a)
{
    if (first_ == a)
    {
        return false;
    }

    DLListBase::link *ap = a->prev_;

    if (ap == first_)
    {
        first_ = a;
        ap->prev_ = a;
    }
    else
    {
        ap->prev_->next_ = a;
    }

    if (a == last_)
    {
        last_ = ap;
        a->next_ = ap;
    }
    else
    {
        a->next_->prev_ = ap;
    }

    a->prev_ = ap->prev_;
    ap->prev_ = a;

    ap->next_ = a->next_;
    a->next_ = ap;

    return true;
}


bool Foam::DLListBase::swapDown(DLListBase::link* a)
{
    if (last_ == a)
    {
        return false;
    }

    DLListBase::link *an = a->next_;

    if (a == first_)
    {
        first_ = an;
        a->prev_ = an;
    }
    else
    {
        a->prev_->next_ = an;
    }

    if (an == last_)
    {
        last_ = a;
        an->next_ = a;
    }
    else
    {
        an->next_->prev_ = a;
    }

    an->prev_ = a->prev_;
    a->prev_ = an;

    a->next_ = an->next_;
    an->next_ = a;

    return true;
}


Foam::DLListBase::link* Foam::DLListBase::removeHead()
{
    if (!first_)
    {
        FatalErrorInFunction
            << "remove from empty list"
            << abort(FatalError);

        // return nullptr;
    }

    --size_;

    DLListBase::link *ret = first_;
    first_ = first_->next_;

    if (!first_)
    {
        last_ = nullptr;
    }

    ret->deregister();
    return ret;
}


Foam::DLListBase::link* Foam::DLListBase::remove(DLListBase::link* item)
{
    --size_;

    DLListBase::link *ret = item;

    if (item == first_ && first_ == last_)
    {
        first_ = nullptr;
        last_ = nullptr;
    }
    else if (item == first_)
    {
        first_ = first_->next_;
        first_->prev_ = first_;
    }
    else if (item == last_)
    {
        last_ = last_->prev_;
        last_->next_ = last_;
    }
    else
    {
        item->next_->prev_ = item->prev_;
        item->prev_->next_ = item->next_;
    }

    ret->deregister();
    return ret;
}


Foam::DLListBase::link* Foam::DLListBase::replace
(
    DLListBase::link* oldLink,
    DLListBase::link* newLink
)
{
    DLListBase::link *ret = oldLink;

    newLink->prev_ = oldLink->prev_;
    newLink->next_ = oldLink->next_;

    if (oldLink == first_ && first_ == last_)
    {
        first_ = newLink;
        last_  = newLink;
    }
    else if (oldLink == first_)
    {
        first_ = newLink;
        newLink->next_->prev_ = newLink;
    }
    else if (oldLink == last_)
    {
        last_ = newLink;
        newLink->prev_->next_ = newLink;
    }
    else
    {
        newLink->prev_->next_ = newLink;
        newLink->next_->prev_ = newLink;
    }

    ret->deregister();
    return ret;
}


// ************************************************************************* //
