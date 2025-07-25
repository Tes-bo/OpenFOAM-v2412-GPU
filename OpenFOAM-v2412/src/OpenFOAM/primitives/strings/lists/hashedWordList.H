/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::hashedWordList

Description
    A wordList with hashed named lookup, which can be faster in some
    situations than using the normal list find/found methods.

SourceFiles
    hashedWordListI.H
    hashedWordList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_hashedWordList_H
#define Foam_hashedWordList_H

#include "wordList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class hashedWordList Declaration
\*---------------------------------------------------------------------------*/

class hashedWordList
:
    public wordList
{
    // Private Data

        //- Lookup HashTable of words vs list-indices
        mutable HashTable<label> lookup_;

public:

    // Constructors

        //- Default construct an empty list
        hashedWordList() = default;

        //- Copy construct
        inline hashedWordList(const hashedWordList& list);

        //- Move construct
        inline hashedWordList(hashedWordList&& list);

        //- Copy construct from list of words
        inline explicit hashedWordList(const wordUList& list);

        //- Copy construct from list of words, eliminating duplicates
        inline hashedWordList(const wordUList& list, bool unique);

        //- Move construct from list of words, optionally eliminating duplicates
        inline hashedWordList(wordList&& list, bool unique=false);

        //- Construct from an initializer list
        inline hashedWordList(std::initializer_list<word> list);

        //- Construct from the word keys of any HashTable, sorting immediately.
        //  This also handles a wordHashSet, which is derived from a HashTable.
        //  The result is similar to a HashTable::sortedToc.
        template<class AnyType, class AnyHash>
        inline explicit hashedWordList
        (
            const HashTable<AnyType, word, AnyHash>& tbl
        );

        //- Construct from Istream
        inline explicit hashedWordList(Istream& is);


    // Member Functions

        //- Clear the list, i.e. set size to zero.
        inline void clear();

        //- Append an element if not already in the list.
        //  \return the change in list length
        inline label push_uniq(const word& val);

        //- Return the hash of words/indices for inspection
        inline const HashTable<label>& lookup() const;

        //- Swap contents
        inline void swap(hashedWordList& list);

        //- Transfer contents of the argument into this list
        //- and annul the argument list, optionally eliminating duplicates
        inline void transfer(hashedWordList& list);

        //- Transfer the contents of the argument List into this list
        //- and annul the argument list, optionally eliminating duplicates
        inline void transfer(wordList& list, bool unique=false);

        //- Rebuild the lookup hash indices
        void rehash() const;

        //- Rebuild the lookup hash indices, or make unique entries first.
        inline void rehash(bool unique);

        //- Adjust the list (if needed) to eliminate duplicate entries,
        //- and rehash the indices
        void uniq();

        //- Inplace sort list and rehash the indices
        inline void sort();


    // Search

        //- Find index of the value (searches the hash).
        //  \return position in list or -1 if not found.
        inline label find(const word& val) const;

        //- Is the value contained in the list (searches the hash).
        inline bool contains(const word& val) const;


    // Member Operators

        //- Return name corresponding to specified index.
        //  Fatal for out of range values.
        inline const word& operator[](const label index) const;

        //- Find index of the value (searches the hash) - same as find().
        //  \return position in list or -1 if not found.
        inline label operator[](const word& val) const;

        //- Check hashed values for the specified name - same as contains().
        //  Can be used as a unary predicate.
        inline bool operator()(const word& val) const;


    // Assignment

        //- Copy assignment. Rehashes the indices.
        inline void operator=(const hashedWordList& list);

        //- Copy assignment from list of words. Rehashes the indices.
        inline void operator=(const wordUList& list);

        //- Copy assignment from initializer list. Rehashes the indices.
        inline void operator=(std::initializer_list<word> list);

        //- Move assignment operator.
        inline void operator=(hashedWordList&& list);

        //- Move assignment from list of words. Rehashes the indices.
        inline void operator=(wordList&& list);


    // Housekeeping

        //- Same as contains(), searches the hash.
        bool found(const word& val) const { return this->contains(val); }

        //- Same as push_uniq()
        FOAM_DEPRECATED_FOR(2022-05, "push_uniq method")
        void append(const word& val) { this->push_uniq(val); }

        //- Same as push_uniq()
        FOAM_DEPRECATED_FOR(2022-10, "push_uniq method")
        void push_back(const word& val) { this->push_uniq(val); }

        //- Same as push_uniq()
        FOAM_DEPRECATED_STRICT(2022-10, "push_uniq method")
        label appendUniq(const word& val) { return this->push_uniq(val); }
};


//- Read from an input stream. Rehashes the indices.
inline Istream& operator>>(Istream& is, hashedWordList& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hashedWordListI.H"

#endif

// ************************************************************************* //
