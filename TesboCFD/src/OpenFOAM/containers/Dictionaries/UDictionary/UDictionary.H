/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::UDictionary

Description
    Template dictionary class which does not manages the storage
    associated with it.

    It is derived from DictionaryBase instantiated on a non-memory managed
    form of intrusive doubly-linked list of \<T\>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_UDictionary_H
#define Foam_UDictionary_H

#include "DictionaryBase.H"
#include "UIDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UDictionary
:
    public DictionaryBase<UIDLList<T>, T>
{
public:

    //- The template instance used for the dictionary content
    typedef DictionaryBase<UIDLList<T>, T> dict_type;


    // Constructors

        //- Default construct: empty without allocation (capacity=0).
        UDictionary() = default;

        //- Construct empty with initial table capacity
        explicit UDictionary(const label initialCapacity)
        :
            dict_type(initialCapacity)
        {}

        //- Copy construct
        UDictionary(const UDictionary& dict)
        :
            dict_type(dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
