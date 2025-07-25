/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "edgeMesh.H"
#include "bitSet.H"
#include "edgeHashes.H"
#include "mergePoints.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(edgeMesh, 0);
    defineRunTimeSelectionTable(edgeMesh, fileExtension);
    defineMemberFunctionSelectionTable(edgeMesh,write,fileExtension);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::wordHashSet Foam::edgeMesh::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


Foam::wordHashSet Foam::edgeMesh::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


bool Foam::edgeMesh::canReadType(const word& fileType, bool verbose)
{
    return checkSupport
    (
        readTypes(),
        fileType,
        verbose,
        "reading"
   );
}


bool Foam::edgeMesh::canWriteType(const word& fileType, bool verbose)
{
    return checkSupport
    (
        writeTypes(),
        fileType,
        verbose,
        "writing"
    );
}


bool Foam::edgeMesh::canRead(const fileName& name, bool verbose)
{
    const word ext =
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );

    return canReadType(ext, verbose);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::edgeMesh::calcPointEdges() const
{
    if (pointEdgesPtr_)
    {
        FatalErrorInFunction
            << "pointEdges already calculated."
            << abort(FatalError);
    }

    pointEdgesPtr_.reset(new labelListList(points_.size()));
    auto& pointEdges = *pointEdgesPtr_;

    invertManyToMany(pointEdges.size(), edges_, pointEdges);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeMesh::clear()
{
    points_.clear();
    edges_.clear();
    pointEdgesPtr_.reset(nullptr);
}


void Foam::edgeMesh::transfer(edgeMesh& mesh)
{
    if (&mesh == this)
    {
        return;  // Self-transfer is a no-op
    }

    points_.transfer(mesh.points_);
    edges_.transfer(mesh.edges_);
    pointEdgesPtr_ = std::move(mesh.pointEdgesPtr_);
}


Foam::label Foam::edgeMesh::regions(labelList& edgeRegion) const
{
    edgeRegion.setSize(edges_.size());
    edgeRegion = -1;

    label startEdgeI = 0;
    label currentRegion = 0;

    while (true)
    {
        while (startEdgeI < edges_.size() && edgeRegion[startEdgeI] != -1)
        {
            startEdgeI++;
        }

        if (startEdgeI == edges_.size())
        {
            break;
        }

        // Found edge that has not yet been assigned a region.
        // Mark connected region with currentRegion starting at startEdgeI.

        edgeRegion[startEdgeI] = currentRegion;
        labelList edgesToVisit(1, startEdgeI);

        while (edgesToVisit.size())
        {
            // neighbours of current edgesToVisit
            DynamicList<label> newEdgesToVisit(edgesToVisit.size());

            // Mark all point connected edges with current region.
            forAll(edgesToVisit, i)
            {
                label edgeI = edgesToVisit[i];

                // Mark connected edges
                const edge& e = edges_[edgeI];

                forAll(e, fp)
                {
                    const labelList& pEdges = pointEdges()[e[fp]];

                    forAll(pEdges, pEdgeI)
                    {
                        label nbrEdgeI = pEdges[pEdgeI];

                        if (edgeRegion[nbrEdgeI] == -1)
                        {
                            edgeRegion[nbrEdgeI] = currentRegion;
                            newEdgesToVisit.append(nbrEdgeI);
                        }
                    }
                }
            }

            edgesToVisit.transfer(newEdgesToVisit);
        }

        currentRegion++;
    }

    return currentRegion;
}


void Foam::edgeMesh::scalePoints(const scalar scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > VSMALL && !equal(scaleFactor, 1))
    {
        points_ *= scaleFactor;
    }
}


void Foam::edgeMesh::mergePoints(const scalar mergeDist)
{
    labelList pointMap;

    label nChanged = Foam::inplaceMergePoints
    (
        points_,
        mergeDist,
        false,
        pointMap
    );

    if (nChanged)
    {
        pointEdgesPtr_.reset(nullptr);   // connectivity change

        for (edge& e : edges_)
        {
            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }

    this->mergeEdges();
}


void Foam::edgeMesh::mergeEdges()
{
    edgeHashSet uniqEdges(2*edges_.size());
    bitSet pointIsUsed(points_.size());

    label nUniqEdges = 0;
    label nUniqPoints = 0;
    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        // Remove degenerate and repeated edges
        // - reordering (e[0] < e[1]) is not really necessary
        if (e[0] != e[1] && uniqEdges.insert(e))
        {
            if (nUniqEdges != edgeI)
            {
                edges_[nUniqEdges] = e;
            }
            edges_[nUniqEdges].sort();
            ++nUniqEdges;

            if (pointIsUsed.set(e[0]))
            {
                ++nUniqPoints;
            }
            if (pointIsUsed.set(e[1]))
            {
                ++nUniqPoints;
            }
        }
    }

    if (debug)
    {
        Info<< "Merging duplicate edges: "
            << (edges_.size() - nUniqEdges)
            << " edges will be deleted, "
            << (points_.size() - nUniqPoints)
            << " unused points will be removed." << endl;
    }

    if (nUniqEdges < edges_.size())
    {
        pointEdgesPtr_.reset(nullptr); // connectivity change
        edges_.setSize(nUniqEdges);  // truncate
    }

    if (nUniqPoints < points_.size())
    {
        pointEdgesPtr_.reset(nullptr); // connectivity change

        // build a oldToNew point-map and rewrite the points.
        // We can do this simultaneously since the point order is unchanged
        // and we are only effectively eliminating some entries.
        labelList pointMap(points_.size(), -1);

        label newId = 0;
        forAll(pointMap, pointi)
        {
            if (pointIsUsed.test(pointi))
            {
                pointMap[pointi] = newId;

                if (newId < pointi)
                {
                    // copy down
                    points_[newId] = points_[pointi];
                }
                ++newId;
            }
        }
        points_.setSize(newId);

        // Renumber edges - already sorted (above)
        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }
}


// ************************************************************************* //
