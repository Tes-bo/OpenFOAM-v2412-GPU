/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::PtrListDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on the memory managed PtrList
    of \<T\> to provide ordered indexing in addition to the dictionary lookup.

SourceFiles
    PtrListDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrListDictionary_H
#define Foam_PtrListDictionary_H

#include "DictionaryBase.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrListDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrListDictionary
:
    public DictionaryBase<PtrList<T>, T>
{
public:

    //- The template instance used for the dictionary content
    typedef DictionaryBase<PtrList<T>, T> dict_type;


    // Constructors

        //- Construct given initial list size
        explicit PtrListDictionary(const label size);

        //- Copy construct
        PtrListDictionary(const PtrListDictionary& dict);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrListDictionary(Istream& is, const INew& inew);

        //- Construct from Istream
        explicit PtrListDictionary(Istream& is);


    // Member functions

        //- Set element to pointer provided and return old element
        autoPtr<T> set(const label i, const word& key, T* ptr);

        //- Set element to autoPtr value provided and return old element
        autoPtr<T> set(const label i, const word& key, autoPtr<T>& ptr);

        //- Set element to tmp value provided and return old element
        autoPtr<T> set(const label i, const word& key, tmp<T>& ptr);


    // Member Operators

        using PtrList<T>::operator[];

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *(dict_type::lookup(key));
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *(dict_type::lookup(key));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrListDictionary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
