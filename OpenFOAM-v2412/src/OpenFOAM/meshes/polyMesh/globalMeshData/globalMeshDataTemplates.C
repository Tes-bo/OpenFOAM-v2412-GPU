/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "globalMeshData.H"
#include "polyMesh.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class CombineOp, class TransformOp>
void Foam::globalMeshData::syncData
(
    List<Type>& elems,
    const labelListList& slaves,
    const labelListList& transformedSlaves,
    const mapDistribute& slavesMap,
    const globalIndexAndTransform& transforms,
    const CombineOp& cop,
    const TransformOp& top
)
{
    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Pull slave data onto master
    slavesMap.distribute(transforms, elems, top);

    // Combine master data with slave data
    forAll(slaves, i)
    {
        Type& elem = elems[i];

        const labelList& slavePoints = slaves[i];

        const labelList& transformSlavePoints =
        (
            transformedSlaves.empty()
          ? labelList::null()
          : transformedSlaves[i]
        );


        // Combine master with untransformed slave data
        for (const label pointi : slavePoints)
        {
            cop(elem, elems[pointi]);
        }

        // Combine master with transformed slave data
        for (const label pointi : transformSlavePoints)
        {
            cop(elem, elems[pointi]);
        }

        // Copy result back to slave slots
        for (const label pointi : slavePoints)
        {
            elems[pointi] = elem;
        }

        for (const label pointi : transformSlavePoints)
        {
            elems[pointi] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute
    (
        transforms,
        elems.size(),
        elems,
        top
    );

    // Reset tag
    UPstream::msgType(oldTag);
}


template<class Type, class CombineOp>
void Foam::globalMeshData::syncData
(
    List<Type>& elems,
    const labelListList& slaves,
    const labelListList& transformedSlaves,
    const mapDistribute& slavesMap,
    const CombineOp& cop
)
{
    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Pull slave data onto master
    slavesMap.distribute(elems);

    // Combine master data with slave data
    forAll(slaves, i)
    {
        Type& elem = elems[i];

        const labelList& slavePoints = slaves[i];

        const labelList& transformSlavePoints =
        (
            transformedSlaves.empty()
          ? labelList::null()
          : transformedSlaves[i]
        );


        // Combine master with untransformed slave data
        for (const label pointi : slavePoints)
        {
            cop(elem, elems[pointi]);
        }

        // Combine master with transformed slave data
        for (const label pointi : transformSlavePoints)
        {
            cop(elem, elems[pointi]);
        }

        // Copy result back to slave slots
        for (const label pointi : slavePoints)
        {
            elems[pointi] = elem;
        }

        for (const label pointi : transformSlavePoints)
        {
            elems[pointi] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems);

    // Reset tag
    UPstream::msgType(oldTag);
}


template<class Type, class CombineOp, class TransformOp>
void Foam::globalMeshData::syncPointData
(
    List<Type>& pointData,
    const CombineOp& cop,
    const TransformOp& top
) const
{
    if (pointData.size() != mesh_.nPoints())
    {
        FatalErrorInFunction
            << "Number of elements in data:" << pointData.size()
            << " differs from number of points in mesh:" << mesh_.nPoints()
            << abort(FatalError);
    }

    // Transfer onto coupled patch
    const indirectPrimitivePatch& cpp = coupledPatch();
    List<Type> cppFld(UIndirectList<Type>(pointData, cpp.meshPoints()));

    syncData
    (
        cppFld,
        globalPointSlaves(),
        globalPointTransformedSlaves(),
        globalPointSlavesMap(),
        globalTransforms(),
        cop,
        top
    );

    // Extract back onto mesh
    forAll(cpp.meshPoints(), i)
    {
        pointData[cpp.meshPoints()[i]] = cppFld[i];
    }
}


// ************************************************************************* //
