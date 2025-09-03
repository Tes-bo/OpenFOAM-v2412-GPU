/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Alias
    Foam::WordPairMap
    A HashTable to objects of type \<T\> with a wordPair key.

Typedef
    Foam::wordPairHashSet

Description
    A HashSet for a wordPair.

Typedef
    Foam::wordPairHashTable

Description
    HashTable of wordPair

\*---------------------------------------------------------------------------*/

#ifndef Foam_wordPairHashes_H
#define Foam_wordPairHashes_H

#include "wordPair.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class T>
    using WordPairMap = HashTable<T, wordPair, Foam::Hash<wordPair>>;

    typedef HashSet<wordPair, Foam::Hash<wordPair>> wordPairHashSet;

    typedef HashTable<word, wordPair, Foam::Hash<wordPair>> wordPairHashTable;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
