/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(const label size)
:
    dict_type(2*size)
{
    PtrList<T>::resize(size);
}


template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(const PtrListDictionary& dict)
:
    dict_type(dict)
{}


template<class T>
template<class INew>
Foam::PtrListDictionary<T>::PtrListDictionary(Istream& is, const INew& iNew)
:
    dict_type(is, iNew)
{}


template<class T>
Foam::PtrListDictionary<T>::PtrListDictionary(Istream& is)
:
    dict_type(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    T* ptr
)
{
    autoPtr<T> old(PtrList<T>::set(i, ptr));

    if (!dict_type::addHashEntry(key, ptr))
    {
        FatalErrorInFunction
            << "Cannot insert with key '" << key << "' into hash-table"
            << abort(FatalError);
    }

    return old;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    autoPtr<T>& ptr
)
{
    return this->set(i, key, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrListDictionary<T>::set
(
    const label i,
    const word& key,
    tmp<T>& ptr
)
{
    return this->set(i, key, ptr.ptr());  // release or clone
}


// ************************************************************************* //
