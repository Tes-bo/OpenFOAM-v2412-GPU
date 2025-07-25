/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::hashedWordList::hashedWordList(const hashedWordList& list)
:
    hashedWordList(static_cast<const wordUList&>(list), false)
{}


inline Foam::hashedWordList::hashedWordList(hashedWordList&& list)
:
    wordList(std::move(static_cast<wordList&>(list))),
    lookup_(std::move(list.lookup_))
{}


inline Foam::hashedWordList::hashedWordList(const wordUList& list)
:
    hashedWordList(list, false)
{}


inline Foam::hashedWordList::hashedWordList(const wordUList& list, bool unique)
:
    wordList(list)
{
    rehash(unique);
}


inline Foam::hashedWordList::hashedWordList(wordList&& list, bool unique)
:
    wordList(std::move(list))
{
    rehash(unique);
}


inline Foam::hashedWordList::hashedWordList(std::initializer_list<word> list)
:
    wordList(list)
{
    rehash();
}


template<class AnyType, class AnyHash>
inline Foam::hashedWordList::hashedWordList
(
    const HashTable<AnyType, word, AnyHash>& tbl
)
:
    wordList(tbl.sortedToc())
{
    rehash();
}


inline Foam::hashedWordList::hashedWordList(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::hashedWordList::clear()
{
    wordList::clear();
    lookup_.clear();
}


inline Foam::label Foam::hashedWordList::push_uniq(const word& val)
{
    if (lookup_.insert(val, size()))
    {
        wordList::push_back(val);
        return 1;  // Increased list length by one
    }

    return 0;  // No change
}


inline const Foam::HashTable<Foam::label>&
Foam::hashedWordList::lookup() const
{
    const label lenList = wordList::size();
    const label lenHash = lookup_.size();

    if ((lenList < lenHash) || (lenList && !lenHash))
    {
        rehash();   // Was somehow out of sync
    }

    return lookup_;
}


inline Foam::label Foam::hashedWordList::find(const word& val) const
{
    return lookup().lookup(val, -1);  // -1 = not found or not hashed
}


inline bool Foam::hashedWordList::contains(const word& val) const
{
    return lookup().contains(val);
}


inline void Foam::hashedWordList::swap(hashedWordList& list)
{
    if (this == &list)
    {
        return;  // Self-swap is a no-op
    }

    wordList::swap(static_cast<wordList&>(list));
    lookup_.swap(list.lookup_);
}


inline void Foam::hashedWordList::transfer(hashedWordList& list)
{
    wordList::transfer(static_cast<wordList&>(list));
    lookup_.transfer(list.lookup_);
}


inline void Foam::hashedWordList::transfer(wordList& list, bool unique)
{
    wordList::transfer(list);
    rehash(unique);
}


inline void Foam::hashedWordList::rehash(bool unique)
{
    if (unique)
    {
        uniq();
    }
    else
    {
        rehash();
    }
}


inline void Foam::hashedWordList::sort()
{
    Foam::sort(*this);
    rehash();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::word& Foam::hashedWordList::operator[]
(
    const label index
) const
{
    return wordList::operator[](index);
}


inline Foam::label Foam::hashedWordList::operator[](const word& val) const
{
    return lookup_.lookup(val, -1);  // -1 = not found or not hashed
}


inline bool Foam::hashedWordList::operator()(const word& val) const
{
    return lookup_.contains(val);
}


inline void Foam::hashedWordList::operator=(const hashedWordList& list)
{
    wordList::operator=(list);
    rehash();
}


inline void Foam::hashedWordList::operator=(const wordUList& list)
{
    wordList::operator=(list);
    rehash();
}


inline void Foam::hashedWordList::operator=(std::initializer_list<word> list)
{
    wordList::operator=(list);
    rehash();
}


inline void Foam::hashedWordList::operator=(hashedWordList&& list)
{
    wordList::transfer(static_cast<wordList&>(list));
    lookup_ = std::move(list.lookup_);
}


inline void Foam::hashedWordList::operator=(wordList&& list)
{
    wordList::transfer(list);
    rehash();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, hashedWordList& list)
{
    is  >> static_cast<wordList&>(list);
    list.rehash();

    return is;
}


// ************************************************************************* //
