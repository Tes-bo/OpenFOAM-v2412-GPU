/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "polyMeshFilter.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class SetType>
void Foam::polyMeshFilter::updateSets(const mapPolyMesh& map)
{
    //
    // Update all sets in memory
    //

    // Note: objectRegistry::lookupClass() instead of
    // objectRegistry::csorted() since it is also used to check
    // contains() in the next bit of code

    const HashTable<const SetType*> sets
    (
        map.mesh().objectRegistry::lookupClass<const SetType>()
    );

    for (const auto& iter : sets.csorted())
    {
        SetType& set = const_cast<SetType&>(*iter.val());
        set.updateMesh(map);
        set.sync(map.mesh());
    }

    //
    // Update all sets on disk
    //

    IOobjectList objs
    (
        map.mesh().time(),
        map.mesh().facesInstance(),
        "polyMesh/sets"
    );

    for (const IOobject& io : objs.csorted<SetType>())
    {
        if (!sets.contains(io.name()))
        {
            // Not in memory. Load it.
            SetType set(io);
            set.updateMesh(map);

            set.write();
        }
    }
}


template<class SetType>
void Foam::polyMeshFilter::copySets
(
    const polyMesh& oldMesh,
    const polyMesh& newMesh
)
{
    for (const SetType& set : oldMesh.objectRegistry::csorted<SetType>())
    {
        auto* setPtr =
            newMesh.objectRegistry::getObjectPtr<SetType>(set.name());

        if (setPtr)
        {
            (*setPtr) = set;
        }
        else
        {
            setPtr = new SetType(newMesh, set.name(), set, set.writeOpt());
            setPtr->store();
        }

        setPtr->sync(newMesh);
    }
}


// ************************************************************************* //
