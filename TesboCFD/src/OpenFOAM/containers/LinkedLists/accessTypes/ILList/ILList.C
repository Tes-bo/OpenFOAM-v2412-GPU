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

#include "ILList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::ILList<LListBase, T>::ILList(const ILList<LListBase, T>& lst)
:
    UILList<LListBase, T>()
{
    for (const auto& item : lst)
    {
        this->push_back(item.clone().ptr());
    }
}


template<class LListBase, class T>
Foam::ILList<LListBase, T>::ILList(ILList<LListBase, T>&& lst)
:
    UILList<LListBase, T>()
{
    LListBase::transfer(lst);
}


template<class LListBase, class T>
template<class CloneArg>
Foam::ILList<LListBase, T>::ILList
(
    const ILList<LListBase, T>& lst,
    const CloneArg& cloneArg
)
:
    UILList<LListBase, T>()
{
    for (const auto& item : lst)
    {
        this->push_back(item.clone(cloneArg).ptr());
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::ILList<LListBase, T>::~ILList()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::ILList<LListBase, T>::pop_front(label n)
{
    if (n > this->size())
    {
        n = this->size();
    }

    while (n > 0)
    {
        T* p = this->removeHead();
        delete p;
        --n;
    }
}


template<class LListBase, class T>
void Foam::ILList<LListBase, T>::clear()
{
    this->pop_front(this->size());
    LListBase::clear();
}


template<class LListBase, class T>
bool Foam::ILList<LListBase, T>::erase(T* item)
{
    T* p = remove(item);
    delete p;
    return bool(p);
}


template<class LListBase, class T>
void Foam::ILList<LListBase, T>::transfer(ILList<LListBase, T>& lst)
{
    clear();
    LListBase::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::ILList<LListBase, T>::operator=(const ILList<LListBase, T>& lst)
{
    this->clear();

    for (const auto& item : lst)
    {
        this->push_back(item.clone().ptr());
    }
}


template<class LListBase, class T>
void Foam::ILList<LListBase, T>::operator=(ILList<LListBase, T>&& lst)
{
    clear();
    LListBase::transfer(lst);
}


// ************************************************************************* //
