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

#include "LList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(const LList<LListBase, T>& lst)
:
    LListBase()
{
    for (const T& val : lst)
    {
        this->push_back(val);
    }
}


template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(LList<LListBase, T>&& lst)
:
    LListBase()
{
    LListBase::transfer(lst);
}


template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(std::initializer_list<T> lst)
:
    LListBase()
{
    for (const T& val : lst)
    {
        this->push_back(val);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LList<LListBase, T>::~LList()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::LList<LListBase, T>::pop_front(label n)
{
    if (n > this->size())
    {
        n = this->size();
    }

    while (n > 0)
    {
        link* p = static_cast<link*>(LListBase::removeHead());
        delete p;
        --n;
    }
}


template<class LListBase, class T>
void Foam::LList<LListBase, T>::clear()
{
    this->pop_front(this->size());
    LListBase::clear();
}


template<class LListBase, class T>
void Foam::LList<LListBase, T>::transfer(LList<LListBase, T>& lst)
{
    clear();
    LListBase::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::LList<LListBase, T>::operator=(const LList<LListBase, T>& lst)
{
    this->clear();

    for (const T& val : lst)
    {
        this->push_back(val);
    }
}


template<class LListBase, class T>
void Foam::LList<LListBase, T>::operator=(LList<LListBase, T>&& lst)
{
    this->clear();

    LListBase::transfer(lst);
}


template<class LListBase, class T>
void Foam::LList<LListBase, T>::operator=(std::initializer_list<T> lst)
{
    this->clear();

    for (const T& val : lst)
    {
        this->push_back(val);
    }
}


// ************************************************************************* //
