/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::EdgeMap

Description
    Map from edge (expressed as its endpoints) to value.
    Hashing (and ==) on an edge is symmetric.

Note
    For easier forward declaration it is implemented as a separate class
    rather than a template alias.

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeHashes_H
#define Foam_edgeHashes_H

#include "edge.H"
#include "HashSet.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Typedefs

//! A HashSet with edge for its key. Hashing (and ==) on an edge is symmetric.
typedef HashSet<edge, Hash<edge>> edgeHashSet;


/*---------------------------------------------------------------------------*\
                           Class EdgeMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class EdgeMap
:
    public HashTable<T, edge, Hash<edge>>
{
public:

    // Constructors

        //- Inherit constructors from HashTable
        using HashTable<T, edge, Hash<edge>>::HashTable;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
