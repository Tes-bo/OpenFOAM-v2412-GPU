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

#include "SLListBase.H"
#include "error.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SLListBase::push_front(SLListBase::link* item)
{
    if (!item)
    {
        return;
    }

    ++size_;

    if (last_)
    {
        item->next_ = last_->next_;
    }
    else
    {
        last_ = item;
    }

    last_->next_ = item;
}


void Foam::SLListBase::push_back(SLListBase::link* item)
{
    if (!item)
    {
        return;
    }

    ++size_;

    if (last_)
    {
        item->next_ = last_->next_;
        last_ = last_->next_ = item;
    }
    else
    {
        last_ = item->next_ = item;
    }
}


Foam::SLListBase::link* Foam::SLListBase::removeHead()
{
    if (last_ == nullptr)
    {
        FatalErrorInFunction
            << "remove from empty list"
            << abort(FatalError);

        // return nullptr;
    }

    --size_;

    SLListBase::link *ret = last_->next_;

    if (ret == last_)
    {
        last_ = nullptr;
    }
    else
    {
        last_->next_ = ret->next_;
    }

    ret->deregister();
    return ret;
}


Foam::SLListBase::link* Foam::SLListBase::remove(SLListBase::link* item)
{
    SLListBase::iterator iter = begin();
    SLListBase::link *prev = iter.get_node();

    if (item == prev)
    {
        return removeHead();
    }

    for (iter.next(); iter != end(); iter.next())
    {
        SLListBase::link *p = iter.get_node();

        if (p == item)
        {
            --size_;

            prev->next_ = p->next_;

            if (p == last_)
            {
                last_ = prev;
            }

            item->deregister();
            return item;
        }

        prev = p;
    }

    // Did not remove
    return nullptr;
}


// ************************************************************************* //
