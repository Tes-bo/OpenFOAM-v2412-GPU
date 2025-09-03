/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "UILList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::UILList<LListBase, T>::UILList(const UILList<LListBase, T>& lst)
{
    for (auto iter = lst.cbegin(); iter != lst.cend(); ++iter)
    {
        this->push_back(&(*iter));
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::UILList<LListBase, T>::operator=(const UILList<LListBase, T>& lst)
{
    LListBase::clear();

    for (auto iter = lst.cbegin(); iter != lst.cend(); ++iter)
    {
        this->push_back(&(*iter));
    }
}


template<class LListBase, class T>
bool Foam::UILList<LListBase, T>::operator==
(
    const UILList<LListBase, T>& rhs
) const
{
    if (this->size() != rhs.size())
    {
        return false;
    }

    auto iter2 = rhs.cbegin();

    for (auto iter1 = this->cbegin(); iter1 != this->cend(); ++iter1, ++iter2)
    {
        if (!(*iter1 == *iter2))
        {
            return false;
        }
    }

    return true;
}


template<class LListBase, class T>
bool Foam::UILList<LListBase, T>::operator!=
(
    const UILList<LListBase, T>& rhs
) const
{
    return !operator==(rhs);
}


// ************************************************************************* //
