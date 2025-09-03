/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::Dictionary

Description
    General purpose template dictionary class that manages the storage
    associated with it.

    It is derived from DictionaryBase instantiated on a memory managed form
    of intrusive doubly-linked list of \<T\>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Dictionary_H
#define Foam_Dictionary_H

#include "DictionaryBase.H"
#include "IDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Dictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Dictionary
:
    public DictionaryBase<IDLList<T>, T>
{
public:

    //- The template instance used for the dictionary content
    typedef DictionaryBase<IDLList<T>, T> dict_type;


    // Constructors

        //- Default construct, or with initial table capacity
        explicit Dictionary(const label initialCapacity = 128)
        :
            dict_type(initialCapacity)
        {}

        //- Copy construct
        Dictionary(const Dictionary& dict)
        :
            dict_type(dict)
        {}


    // Member Functions

        //- Remove an entry specified by keyword and delete the pointer.
        //  \return true if the keyword was found
        bool erase(const word& keyword)
        {
            T* ptr = this->remove(keyword);
            delete ptr;
            return bool(ptr);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
