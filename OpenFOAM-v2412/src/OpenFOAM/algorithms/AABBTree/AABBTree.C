/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "AABBTree.H"
#include "bitSet.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::AABBTree<Type>::writeOBJ
(
    const bool leavesOnly,
    const bool writeLinesOnly,
    const treeBoundBox& bb,
    const label nodeI,
    const List<Pair<treeBoundBox>>& bbs,
    const List<Pair<label>>& nodes,
    label& vertI,
    Ostream& os
) const
{
    if (!leavesOnly || nodeI < 0)
    {
        AABBTreeBase::writeOBJ(os, bb, vertI, writeLinesOnly);
    }

    // recurse to find leaves
    if (nodeI >= 0)
    {
        writeOBJ
        (
            leavesOnly,
            writeLinesOnly,
            bbs[nodeI].first(),
            nodes[nodeI].first(),
            bbs,
            nodes,
            vertI,
            os
        );
        writeOBJ
        (
            leavesOnly,
            writeLinesOnly,
            bbs[nodeI].second(),
            nodes[nodeI].second(),
            bbs,
            nodes,
            vertI,
            os
        );
    }
}


template<class Type>
void Foam::AABBTree<Type>::createBoxes
(
    const bool equalBinSize,
    const label level,
    const UList<Type>& objects,
    const pointField& points,
    const labelUList& objectIDs,
    const treeBoundBox& bb,
    const label nodeI,

    DynamicList<Pair<treeBoundBox>>& bbs,
    DynamicList<labelPair>& nodes,
    DynamicList<labelList>& addressing
) const
{
    // Determine which direction to divide the box

    const direction maxDir = bb.maxDir();
    const scalar maxSpan = bb.span()[maxDir];

    scalar pivotValue;

    if (equalBinSize)
    {
        // Pick up points used by this set of objects

        bitSet isUsedPoint(points.size());
        DynamicList<scalar> component(points.size());

        for (const label objI : objectIDs)
        {
            const Type& obj = objects[objI];

            for (const label pointI : obj)
            {
                if (isUsedPoint.set(pointI))
                {
                    component.push_back(points[pointI][maxDir]);
                }
            }
        }

        // Determine the median

        Foam::sort(component);

        pivotValue = component[component.size()/2];
    }
    else
    {
        // Geometric middle
        pivotValue = bb.min()[maxDir] + 0.5*maxSpan;
    }


    const scalar divMin = pivotValue + tolerance_*maxSpan;
    const scalar divMax = pivotValue - tolerance_*maxSpan;


    // Assign the objects to min or max bin

    DynamicList<label> minBinObjectIDs(objectIDs.size());
    DynamicList<label> maxBinObjectIDs(objectIDs.size());

    treeBoundBox minBb;
    treeBoundBox maxBb;

    for (const label objI : objectIDs)
    {
        const Type& obj = objects[objI];

        bool addMin = false;
        bool addMax = false;

        for (const label pointi : obj)
        {
            const scalar& cmptValue = points[pointi][maxDir];

            addMin = addMin || (cmptValue < divMin);
            addMax = addMax || (cmptValue > divMax);

            if (addMin && addMax) break;
        }

        // Note: object is inserted into both min/max child boxes (duplicated)
        // if it crosses the bin boundaries
        if (addMin)
        {
            minBinObjectIDs.push_back(objI);
            minBb.add(points, obj);
        }
        if (addMax)
        {
            maxBinObjectIDs.push_back(objI);
            maxBb.add(points, obj);
        }
    }

    // Inflate box in case geometry reduces to 2-D
    if (minBinObjectIDs.size())
    {
        minBb.inflate(0.01);
    }
    if (maxBinObjectIDs.size())
    {
        maxBb.inflate(0.01);
    }

    minBinObjectIDs.shrink();
    maxBinObjectIDs.shrink();

    bool addMin = (minBinObjectIDs.size() > minLeafSize_ && level < maxLevel_);
    bool addMax = (maxBinObjectIDs.size() > minLeafSize_ && level < maxLevel_);

    // Since bounding boxes overlap, verify that splitting was effective

    if
    (
        objectIDs.size() <= (minBinObjectIDs.size() + minLeafSize_/2)
     || objectIDs.size() <= (maxBinObjectIDs.size() + minLeafSize_/2)
    )
    {
        addMin = addMax = false;
    }

    label minI;
    if (addMin)
    {
        // New leaf
        minI = nodes.size();
        nodes.emplace_back(-1, -1);
    }
    else
    {
        // Update existing leaf
        minI = -addressing.size() - 1;
        addressing.push_back(minBinObjectIDs);
    }

    label maxI;
    if (addMax)
    {
        // New leaf
        maxI = nodes.size();
        nodes.emplace_back(-1, -1);
    }
    else
    {
        // Update existing leaf
        maxI = -addressing.size() - 1;
        addressing.push_back(maxBinObjectIDs);
    }

    nodes(nodeI) = labelPair(minI, maxI);
    bbs(nodeI) = Pair<treeBoundBox>(minBb, maxBb);

    // Recurse
    if (minI >= 0)
    {
        createBoxes
        (
            equalBinSize,
            level + 1,
            objects,
            points,
            minBinObjectIDs,
            minBb,
            minI,
            bbs,
            nodes,
            addressing
        );
    }
    if (maxI >= 0)
    {
        createBoxes
        (
            equalBinSize,
            level + 1,
            objects,
            points,
            maxBinObjectIDs,
            maxBb,
            maxI,
            bbs,
            nodes,
            addressing
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AABBTree<Type>::AABBTree()
:
    maxLevel_(0),
    minLeafSize_(0),
    boundBoxes_(),
    addressing_()
{}


template<class Type>
Foam::AABBTree<Type>::AABBTree
(
    const UList<Type>& objects,
    const pointField& points,
    const bool equalBinSize,
    label maxLevel,
    label minLeafSize
)
:
    maxLevel_(maxLevel),
    minLeafSize_(minLeafSize),
    boundBoxes_(),
    addressing_()
{
    if (objects.empty())
    {
        return;
    }


    DynamicList<Pair<treeBoundBox>> bbs(maxLevel);
    DynamicList<labelPair> nodes(maxLevel);
    DynamicList<labelList> addr(maxLevel);

    nodes.emplace_back(-1, -1);
    treeBoundBox topBb(points);
    topBb.inflate(0.01);

    labelList objectIDs(identity(objects.size()));

    createBoxes
    (
        equalBinSize,
        0,          // starting at top level
        objects,
        points,
        objectIDs,
        topBb,
        0,          // starting node

        bbs,
        nodes,
        addr
    );


    //{
    //    OFstream os("tree.obj");
    //    label vertI = 0;
    //    writeOBJ
    //    (
    //        true,       // leavesOnly
    //        false,      // writeLinesOnly
    //
    //        topBb,
    //        0,
    //        bbs,
    //        nodes,
    //        vertI,
    //        os
    //    );
    //}


    // transfer flattened tree to persistent storage
    DynamicList<treeBoundBox> boundBoxes(2*bbs.size());
    DynamicList<labelList> addressing(2*addr.size());

    forAll(nodes, nodeI)
    {
        if (nodes[nodeI].first() < 0)
        {
            boundBoxes.push_back(bbs[nodeI].first());
            addressing.push_back(addr[-(nodes[nodeI].first() + 1)]);
        }
        if (nodes[nodeI].second() < 0)
        {
            boundBoxes.push_back(bbs[nodeI].second());
            addressing.push_back(addr[-(nodes[nodeI].second() + 1)]);
        }
    }

    boundBoxes_.transfer(boundBoxes);
    addressing_.transfer(addressing);


    if (0)
    {
        bitSet checked(objects.size());
        for (const auto& box : addressing_)
        {
            for (const auto& id : box)
            {
                checked.set(id);
            }
        }

        const label unsetSize = checked.count(false);

        if (unsetSize)
        {
            Info<< "*** Problem: IDs not set: " << unsetSize << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AABBTree<Type>::writeOBJ(Ostream& os) const
{
    label vertIndex(0);

    for (const treeBoundBox& bb : boundBoxes_)
    {
        // writeLinesOnly=false
        AABBTreeBase::writeOBJ(os, bb, vertIndex, false);
    }
}


template<class Type>
bool Foam::AABBTree<Type>::pointInside(const point& pt) const
{
    for (const treeBoundBox& bb : boundBoxes_)
    {
        if (bb.contains(pt))
        {
            return true;
        }
    }

    return false;
}


template<class Type>
bool Foam::AABBTree<Type>::overlaps(const boundBox& bbIn) const
{
    for (const treeBoundBox& bb : boundBoxes_)
    {
        if (bb.overlaps(bbIn))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const AABBTree<Type>& tree)
{
    os  << tree.boundBoxes_ << tree.addressing_;

    os.check(FUNCTION_NAME);
    return os;
}


template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, AABBTree<Type>& tree)
{
    is  >> tree.boundBoxes_ >> tree.addressing_;

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
