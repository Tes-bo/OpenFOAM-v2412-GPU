/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PtrDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on a memory managed form of
    intrusive doubly-linked list of \<T\>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrDictionary_H
#define Foam_PtrDictionary_H

#include "DictionaryBase.H"
#include "DLPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrDictionary
:
    public DictionaryBase<DLPtrList<T>, T>
{
public:

    //- The template instance used for the dictionary content
    typedef DictionaryBase<DLPtrList<T>, T> dict_type;


    // Constructors

        //- Default construct: empty without allocation (capacity=0).
        PtrDictionary() = default;

        //- Construct empty with initial table capacity
        explicit PtrDictionary(const label initialCapacity)
        :
            dict_type(initialCapacity)
        {}

        //- Copy construct
        PtrDictionary(const PtrDictionary& dict)
        :
            dict_type(dict)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrDictionary(Istream& is, const INew& inew)
        :
            dict_type(is, inew)
        {}

        //- Construct from Istream
        explicit PtrDictionary(Istream& is)
        :
            dict_type(is)
        {}


    // Member Operators

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

#endif

// ************************************************************************* //
