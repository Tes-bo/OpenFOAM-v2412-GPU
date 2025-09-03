/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::DynamicID

Description
    A class that holds the data needed to identify things (zones, patches)
    in a dynamic mesh.

    The thing is identified by name.
    Its indices are updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef DynamicID_H
#define DynamicID_H

#include "wordRe.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class DynamicID Declaration
\*---------------------------------------------------------------------------*/

template<class ObjectType>
class DynamicID
{
    // Private Data

        //- Selector name
        wordRe key_;

        //- Selection indices
        labelList indices_;


public:

    // Constructors

        //- Construct from selector name and object
        DynamicID(const wordRe& key, const ObjectType& obj)
        :
            key_(key),
            indices_(obj.indices(key_))
        {}

        //- Construct from selector name and object
        DynamicID(wordRe&& key, const ObjectType& obj)
        :
            key_(std::move(key)),
            indices_(obj.indices(key_))
        {}

        //- Construct from selector name and object
        DynamicID(const word& key, const ObjectType& obj)
        :
            DynamicID(wordRe(key), obj)
        {}

        //- Construct from selector name and object
        DynamicID(const keyType& key, const ObjectType& obj)
        :
            DynamicID(wordRe(key), obj)
        {}

        //- Construct from Istream and object
        DynamicID(Istream& is, const ObjectType& obj)
        :
            DynamicID(wordRe(is), obj)
        {}


    //- Destructor
    ~DynamicID() = default;


    // Member Functions

    // Access

        //- The selector name
        const wordRe& name() const noexcept
        {
            return key_;
        }

        //- The indices of matching items
        const labelList& indices() const noexcept
        {
            return indices_;
        }

        //- The index of the first matching items, -1 if no matches
        label index() const
        {
            return indices_.empty() ? -1 : indices_.first();
        }

        //- Has the zone been found
        bool active() const noexcept
        {
            return !indices_.empty();
        }


    // Edit

        //- Update
        void update(const ObjectType& obj)
        {
            indices_ = obj.indices(key_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ObjectType>
Ostream& operator<<(Ostream& os, const DynamicID<ObjectType>& obj)
{
    os  << token::BEGIN_LIST
        << obj.name() << token::SPACE << obj.index()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
