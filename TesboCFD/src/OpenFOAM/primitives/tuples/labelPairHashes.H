/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::LabelPairMap

Description
    A HashTable to objects of type \<T\> with a labelPair key.
    The hashing is based on labelPair (FixedList) and is thus non-commutative.
    For commutative hashing, an EdgeMap can be used.

Typedef
    Foam::labelPairHashSet

Description
    A HashSet for a labelPair.
    The hashing is based on labelPair (FixedList) and is thus non-commutative.

Typedef
    Foam::labelPairLookup

Description
    This is a Map of a labelPair to a label.
    Used for e.g. for face1, face2 to shared edge.
    The hashing is based on labelPair (FixedList) and is thus non-commutative.

\*---------------------------------------------------------------------------*/

#ifndef labelPairHashes_H
#define labelPairHashes_H

#include "labelPair.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class T>
    using LabelPairMap = HashTable<T, labelPair, Foam::Hash<labelPair>>;

    typedef HashSet<labelPair, Foam::Hash<labelPair>> labelPairHashSet;

    typedef HashTable<label, labelPair, Foam::Hash<labelPair>> labelPairLookup;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
