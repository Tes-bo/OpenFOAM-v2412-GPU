/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
   This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Subsetting of an adjacency matrix (as CompactListList).
    Can be relocated elsewhere.

\*---------------------------------------------------------------------------*/

#include "CompactListList.H"
#include "bitSet.H"
#include "ListOps.H"
#include "Map.H"

namespace Foam
{

// Perform a subset of the adjacency matrix
CompactListList<label> subsetAdjacency
(
    const bitSet& select,   // could also be labelHashSet
    const CompactListList<label>& input,
    labelList& subMap
)
{
    // Corresponds to cellMap etc (the original selection)
    subMap = select.sortedToc();

    // Ensure that the subMap corresponds to a valid subset
    {
        label validSize = 0;

        const label nTotal = input.size();

        forAllReverse(subMap, i)
        {
            if (subMap[i] < nTotal)
            {
                validSize = i + 1;
                break;
            }
        }

        subMap.resize(validSize);
    }


    // Assumed to be sparse - use Map for reverse lookup
    const Map<label> reverseMap(invertToMap(subMap));


    // Pass 1: determine the selected sub-sizes
    labelList sizes(subMap.size(), Foam::zero{});

    forAll(subMap, idx)
    {
        for (const label nbr : input[subMap[idx]])
        {
            if
            (
                select.test(nbr)
             && reverseMap.contains(nbr)  // extra consistency (paranoid)
            )
            {
                ++sizes[idx];
            }
        }
    }


    CompactListList<label> output(sizes);

    // Reuse sizes as output offset into output.values()
    sizes = labelList::subList(output.offsets(), output.size());
    labelList& values = output.values();


    // Pass 2: extract sub-adjacent matrix

    label newNbr = -1;

    forAll(subMap, idx)
    {
        for (const label nbr : input[subMap[idx]])
        {
            if
            (
                select.test(nbr)
             && (newNbr = reverseMap.lookup(nbr, -1)) >= 0
            )
            {
                values[sizes[idx]++] = newNbr;
            }
        }
    }

    return output;
}


// Perform a subset of the adjacency matrix
CompactListList<label> subsetAdjacency
(
    const labelRange& slice,
    const CompactListList<label>& input,
    labelList& subMap
)
{
    // Ensure that the selection corresponds to a valid subset
    const labelRange select = slice.subset0(input.size());

    // Corresponds to cellMap etc (the original selection)
    subMap = Foam::identity(select);


    // Pass 1: determine the selected sub-sizes
    labelList sizes(subMap.size(), Foam::zero{});

    forAll(subMap, idx)
    {
        for (const label nbr : input[subMap[idx]])
        {
            if (select.contains(nbr))
            {
                ++sizes[idx];
            }
        }
    }


    CompactListList<label> output(sizes);

    // Reuse sizes as output offset into output.values()
    sizes = labelList::subList(output.offsets(), output.size());
    labelList& values = output.values();


    // Pass 2: extract sub-adjacent matrix

    const label localOffset = select.start();

    forAll(subMap, idx)
    {
        for (const label nbr : input[subMap[idx]])
        {
            if (select.contains(nbr))
            {
                values[sizes[idx]++] = nbr - localOffset;
            }
        }
    }

    return output;
}

} // End namespace Foam


// ************************************************************************* //
