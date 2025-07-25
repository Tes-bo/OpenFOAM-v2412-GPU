/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::topoBitSet

Description
    Base for a special purpose topoSet using labels stored as a bitSet.

SourceFiles
    topoBitSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_topoBitSet_H
#define Foam_topoBitSet_H

#include "topoSet.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class topoBitSet Declaration
\*---------------------------------------------------------------------------*/

class topoBitSet
:
    public topoSet
{
protected:

    // Protected Data

        bitSet selected_;


    // Protected Member Functions

        //- Update map from map.
        //  Used to update cell/face labels after morphing
        virtual void updateLabels(const labelUList& map);

        //- Check limits on addressable range.
        virtual void check(const label maxSize);


        //- Construct (no-read) with empty selection
        topoBitSet(const polyMesh& mesh, const word& setName);

        //- Construct (no-read) with \c size elements of initial value
        topoBitSet
        (
            const polyMesh& mesh,
            const word& setName,
            const label size,
            const bool val
        );

        //- Copy construct from bitset values, size elements
        topoBitSet
        (
            const polyMesh& mesh,
            const word& setName,
            const label size,
            const bitSet& bits
        );

        //- Move construct from bitset values, size elements
        topoBitSet
        (
            const polyMesh& mesh,
            const word& setName,
            const label size,
            bitSet&& bits
        );


public:

    //- Destructor
    virtual ~topoBitSet() = default;


    // Member Functions

        //- Return the bitSet
        const bitSet& addressing() const noexcept
        {
            return selected_;
        }

        //- Access the bitSet
        bitSet& addressing() noexcept
        {
            return selected_;
        }

        //- Set values to false, leaving the size untouched
        void reset()
        {
            selected_.reset();
        }

        //- Has the given index?
        virtual bool contains(const label id) const;

        //- Has the given index?
        virtual bool found(const label id) const;

        //- Set an index
        virtual bool set(const label id);

        //- Unset an index
        virtual bool unset(const label id);

        //- Set multiple indices
        virtual void set(const labelUList& labels);

        //- Unset multiple indices
        virtual void unset(const labelUList& labels);

        //- Invert contents.
        //  Insert all members [0,maxLen) which were not in set.
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const labelUList& elems);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add given elements to the set
        virtual void addSet(const labelUList& elems);

        //- Add given elements to the set
        virtual void addSet(const topoSet& set);

        //- Subtract given elements from the set
        virtual void subtractSet(const labelUList& elems);

        //- Subtract given elements from the set
        virtual void subtractSet(const topoSet& set);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
