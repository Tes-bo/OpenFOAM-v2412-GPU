/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "extendedEdgeMesh.H"
#include "surfaceFeatures.H"
#include "triSurface.H"
#include "Random.H"
#include "Time.H"
#include "OBJstream.H"
#include "DynamicField.H"
#include "edgeMeshFormatsCore.H"
#include "IOmanip.H"
#include "searchableSurface.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedEdgeMesh, 0);
}


const Foam::Enum
<
    Foam::extendedEdgeMesh::pointStatus
>
Foam::extendedEdgeMesh::pointStatusNames_
({
    { pointStatus::CONVEX, "convex" },
    { pointStatus::CONCAVE, "concave" },
    { pointStatus::MIXED, "mixed" },
    { pointStatus::NONFEATURE, "nonFeature" },
});


const Foam::Enum
<
    Foam::extendedEdgeMesh::edgeStatus
>
Foam::extendedEdgeMesh::edgeStatusNames_
({
    { edgeStatus::EXTERNAL, "external" },
    { edgeStatus::INTERNAL, "internal" },
    { edgeStatus::FLAT, "flat" },
    { edgeStatus::OPEN, "open" },
    { edgeStatus::MULTIPLE, "multiple" },
    { edgeStatus::NONE, "none" },
});


const Foam::Enum
<
    Foam::extendedEdgeMesh::sideVolumeType
>
Foam::extendedEdgeMesh::sideVolumeTypeNames_
({
    { sideVolumeType::INSIDE, "inside" },
    { sideVolumeType::OUTSIDE, "outside" },
    { sideVolumeType::BOTH, "both" },
    { sideVolumeType::NEITHER, "neither" },
});


Foam::scalar Foam::extendedEdgeMesh::cosNormalAngleTol_ =
    Foam::cos(degToRad(0.1));


Foam::label Foam::extendedEdgeMesh::convexStart_ = 0;

Foam::label Foam::extendedEdgeMesh::externalStart_ = 0;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::wordHashSet Foam::extendedEdgeMesh::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


Foam::wordHashSet Foam::extendedEdgeMesh::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


bool Foam::extendedEdgeMesh::canReadType(const word& fileType, bool verbose)
{
    return edgeMeshFormatsCore::checkSupport
    (
        readTypes(),
        fileType,
        verbose,
        "reading"
   );
}


bool Foam::extendedEdgeMesh::canWriteType(const word& fileType, bool verbose)
{
    return edgeMeshFormatsCore::checkSupport
    (
        writeTypes(),
        fileType,
        verbose,
        "writing"
    );
}


bool Foam::extendedEdgeMesh::canRead(const fileName& name, bool verbose)
{
    const word ext =
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );

    return canReadType(ext, verbose);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::extendedEdgeMesh::pointStatus
Foam::extendedEdgeMesh::classifyFeaturePoint
(
    label ptI
) const
{
    const labelList& ptEds(pointEdges()[ptI]);

    label nPtEds = ptEds.size();
    label nExternal = 0;
    label nInternal = 0;

    if (nPtEds == 0)
    {
        // There are no edges attached to the point, this is a problem
        return NONFEATURE;
    }

    forAll(ptEds, i)
    {
        edgeStatus edStat = getEdgeStatus(ptEds[i]);

        if (edStat == EXTERNAL)
        {
            nExternal++;
        }
        else if (edStat == INTERNAL)
        {
            nInternal++;
        }
    }

    if (nExternal == nPtEds)
    {
        return CONVEX;
    }
    else if (nInternal == nPtEds)
    {
        return CONCAVE;
    }

    return MIXED;
}


void Foam::extendedEdgeMesh::cut
(
    const searchableSurface& surf,

    labelList& pointMap,
    labelList& edgeMap,
    labelList& pointsFromEdge,
    labelList& oldEdge,
    labelList& surfTri
)
{
    const edgeList& edges = this->edges();
    const pointField& points = this->points();


    List<List<pointIndexHit>> edgeHits(edges.size());
    {
        pointField start(edges.size());
        pointField end(edges.size());
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];
            start[edgeI] = points[e[0]];
            end[edgeI] = points[e[1]];
        }
        surf.findLineAll(start, end, edgeHits);
    }

    // Count number of hits
    label nHits = 0;
    forAll(edgeHits, edgeI)
    {
        nHits += edgeHits[edgeI].size();
    }

    DynamicField<point> newPoints(points);
    DynamicList<label> newToOldPoint(identity(points.size()));

    newPoints.setCapacity(newPoints.size()+nHits);
    newToOldPoint.setCapacity(newPoints.capacity());

    DynamicList<edge> newEdges(edges);
    DynamicList<label> newToOldEdge(identity(edges.size()));

    newEdges.setCapacity(newEdges.size()+nHits);
    newToOldEdge.setCapacity(newEdges.capacity());

    // Information on additional points
    DynamicList<label> dynPointsFromEdge(nHits);
    DynamicList<label> dynOldEdge(nHits);
    DynamicList<label> dynSurfTri(nHits);

    forAll(edgeHits, edgeI)
    {
        const List<pointIndexHit>& eHits = edgeHits[edgeI];

        if (eHits.size())
        {
            label prevPtI = edges[edgeI][0];
            forAll(eHits, eHitI)
            {
                label newPtI = newPoints.size();

                newPoints.append(eHits[eHitI].hitPoint());
                newToOldPoint.append(edges[edgeI][0]);  // map from start point
                dynPointsFromEdge.append(newPtI);
                dynOldEdge.append(edgeI);
                dynSurfTri.append(eHits[eHitI].index());

                if (eHitI == 0)
                {
                    newEdges[edgeI] = edge(prevPtI, newPtI);
                }
                else
                {
                    newEdges.append(edge(prevPtI, newPtI));
                    newToOldEdge.append(edgeI);
                }
                prevPtI = newPtI;
            }
            newEdges.append(edge(prevPtI, edges[edgeI][1]));
            newToOldEdge.append(edgeI);
        }
    }

    pointField allPoints;
    allPoints.transfer(newPoints);
    pointMap.transfer(newToOldPoint);

    edgeList allEdges;
    allEdges.transfer(newEdges);
    edgeMap.transfer(newToOldEdge);

    pointsFromEdge.transfer(dynPointsFromEdge);
    oldEdge.transfer(dynOldEdge);
    surfTri.transfer(dynSurfTri);

    // Update local information
    autoMap(allPoints, allEdges, pointMap, edgeMap);
}


void Foam::extendedEdgeMesh::select
(
    const searchableSurface& surf,
    const volumeType volType,   // side to keep
    labelList& pointMap,        // sub to old points
    labelList& edgeMap          // sub to old edges
)
{
    const edgeList& edges = this->edges();
    const pointField& points = this->points();

    // Test edge centres for inside/outside
    if (volType == volumeType::INSIDE || volType == volumeType::OUTSIDE)
    {
        pointField edgeCentres(edges.size());
        forAll(edgeCentres, edgeI)
        {
            const edge& e = edges[edgeI];
            edgeCentres[edgeI] = e.centre(points);
        }
        List<volumeType> volTypes;
        surf.getVolumeType(edgeCentres, volTypes);

        // Extract edges on correct side
        edgeMap.setSize(edges.size());
        label compactEdgeI = 0;

        forAll(volTypes, edgeI)
        {
            if (volTypes[edgeI] == volType)
            {
                edgeMap[compactEdgeI++] = edgeI;
            }
        }
        edgeMap.setSize(compactEdgeI);

        // Extract used points
        labelList pointToCompact(points.size(), -1);
        forAll(edgeMap, i)
        {
            const edge& e = edges[edgeMap[i]];
            pointToCompact[e[0]] = labelMax;       // tag with a value
            pointToCompact[e[1]] = labelMax;
        }

        pointMap.setSize(points.size());
        label compactPointI = 0;
        forAll(pointToCompact, pointI)
        {
            if (pointToCompact[pointI] != -1)
            {
                pointToCompact[pointI] = compactPointI;
                pointMap[compactPointI++] = pointI;
            }
        }
        pointMap.setSize(compactPointI);
        pointField subPoints(points, pointMap);

        // Renumber edges
        edgeList subEdges(edgeMap.size());
        forAll(edgeMap, i)
        {
            const edge& e = edges[edgeMap[i]];
            subEdges[i][0] = pointToCompact[e[0]];
            subEdges[i][1] = pointToCompact[e[1]];
        }

        // Reset primitives and map other quantities
        autoMap(subPoints, subEdges, pointMap, edgeMap);
    }
    else
    {
        pointMap = identity(points.size());
        edgeMap = identity(edges.size());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedEdgeMesh::extendedEdgeMesh()
:
    edgeMesh(),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(nullptr),
    edgeTree_(nullptr),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh(std::nullptr_t)
:
    edgeMesh(),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(nullptr),
    edgeTree_(nullptr),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh(const extendedEdgeMesh& fem)
:
    edgeMesh(fem),
    concaveStart_(fem.concaveStart()),
    mixedStart_(fem.mixedStart()),
    nonFeatureStart_(fem.nonFeatureStart()),
    internalStart_(fem.internalStart()),
    flatStart_(fem.flatStart()),
    openStart_(fem.openStart()),
    multipleStart_(fem.multipleStart()),
    normals_(fem.normals()),
    normalVolumeTypes_(fem.normalVolumeTypes()),
    edgeDirections_(fem.edgeDirections()),
    normalDirections_(fem.normalDirections()),
    edgeNormals_(fem.edgeNormals()),
    featurePointNormals_(fem.featurePointNormals()),
    featurePointEdges_(fem.featurePointEdges()),
    regionEdges_(fem.regionEdges()),
    pointTree_(nullptr),
    edgeTree_(nullptr),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh(Istream& is)
{
    is >> *this;
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const pointField& points,
    const edgeList& edges
)
:
    extendedEdgeMesh()
{
    this->storedPoints() = points;
    this->storedEdges()  = edges;
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    pointField&& points,
    edgeList&& edges
)
:
    extendedEdgeMesh()
{
    this->storedPoints().transfer(points);
    this->storedEdges().transfer(edges);
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const surfaceFeatures& sFeat,
    const boolList& surfBaffleRegions
)
:
    extendedEdgeMesh(nullptr)
{
    // Extract and reorder the data from surfaceFeatures
    const triSurface& surf = sFeat.surface();
    const labelList& featureEdges = sFeat.featureEdges();
    const labelList& featurePoints = sFeat.featurePoints();

    // Get a labelList of all the featureEdges that are region edges
    const labelList regionFeatureEdges(identity(sFeat.nRegionEdges()));

    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );

    const labelListList& edgeFaces = surf.edgeFaces();

    normalVolumeTypes_.setSize(normals_.size());

    // Noting when the normal of a face has been used so not to duplicate
    labelList faceMap(surf.size(), -1);

    label nAdded = 0;

    forAll(featureEdges, i)
    {
        label sFEI = featureEdges[i];

        // Pick up the faces adjacent to the feature edge
        const labelList& eFaces = edgeFaces[sFEI];

        forAll(eFaces, j)
        {
            label eFI = eFaces[j];

            // Check to see if the points have been already used
            if (faceMap[eFI] == -1)
            {
                normalVolumeTypes_[nAdded++] =
                    (
                        surfBaffleRegions[surf[eFI].region()]
                      ? BOTH
                      : INSIDE
                    );

                faceMap[eFI] = nAdded - 1;
            }
        }
    }
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const PrimitivePatch<faceList, pointField>& surf,
    const labelUList& featureEdges,
    const labelUList& regionFeatureEdges,
    const labelUList& featurePoints
)
:
    extendedEdgeMesh(nullptr)
{
    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const pointField& pts,
    const edgeList& eds,
    label concaveStart,
    label mixedStart,
    label nonFeatureStart,
    label internalStart,
    label flatStart,
    label openStart,
    label multipleStart,
    const vectorField& normals,
    const List<sideVolumeType>& normalVolumeTypes,
    const vectorField& edgeDirections,
    const labelListList& normalDirections,
    const labelListList& edgeNormals,
    const labelListList& featurePointNormals,
    const labelListList& featurePointEdges,
    const labelList& regionEdges
)
:
    edgeMesh(pts, eds),
    concaveStart_(concaveStart),
    mixedStart_(mixedStart),
    nonFeatureStart_(nonFeatureStart),
    internalStart_(internalStart),
    flatStart_(flatStart),
    openStart_(openStart),
    multipleStart_(multipleStart),
    normals_(normals),
    normalVolumeTypes_(normalVolumeTypes),
    edgeDirections_(edgeDirections),
    normalDirections_(normalDirections),
    edgeNormals_(edgeNormals),
    featurePointNormals_(featurePointNormals),
    featurePointEdges_(featurePointEdges),
    regionEdges_(regionEdges),
    pointTree_(nullptr),
    edgeTree_(nullptr),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const fileName& name,
    const word& fileType
)
:
    extendedEdgeMesh()
{
    read(name, fileType);
}


Foam::extendedEdgeMesh::extendedEdgeMesh(const fileName& name)
:
    extendedEdgeMesh()
{
    read(name);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::extendedEdgeMesh::read(const fileName& name)
{
    if (name.has_ext("gz"))
    {
        return read(name.lessExt(), name.stem().ext());
    }

    return read(name, name.ext());
}


bool Foam::extendedEdgeMesh::read
(
    const fileName& name,
    const word& fileType
)
{
    // Read via selector mechanism
    transfer(*New(name, fileType));
    return true;
}


void Foam::extendedEdgeMesh::nearestFeaturePoint
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = pointTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedEdgeMesh::nearestFeatureEdge
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = edgeTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedEdgeMesh::nearestFeatureEdge
(
    const pointField& samples,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        nearestFeatureEdge
        (
            samples[i],
            searchDistSqr[i],
            info[i]
        );
    }
}


void Foam::extendedEdgeMesh::nearestFeatureEdgeByType
(
    const point& sample,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    const PtrList<indexedOctree<treeDataEdge>>& edgeTrees = edgeTreesByType();

    info.resize(edgeTrees.size());

    forAll(edgeTrees, i)
    {
        const auto& tree = edgeTrees[i];
        const auto& treeData = edgeTrees[i].shapes();

        info[i] = tree.findNearest(sample, searchDistSqr[i]);

        // The index returned by the indexedOctree is local to the slice of
        // edges it was supplied with, return the index to the value in the
        // complete edge list

        info[i].setIndex(treeData.objectIndex(info[i].index()));
    }
}


void Foam::extendedEdgeMesh::allNearestFeaturePoints
(
    const point& sample,
    scalar searchRadiusSqr,
    List<pointIndexHit>& info
) const
{
    // Pick up all the feature points that intersect the search sphere
    labelList elems = pointTree().findSphere
    (
        sample,
        searchRadiusSqr
    );

    DynamicList<pointIndexHit> dynPointHit(elems.size());

    const auto& treeData = pointTree().shapes();

    for (const label index : elems)
    {
        const point& pt = treeData.centre(index);

        pointIndexHit nearHit(true, pt, index);

        dynPointHit.append(nearHit);
    }

    info.transfer(dynPointHit);
}


void Foam::extendedEdgeMesh::allNearestFeatureEdges
(
    const point& sample,
    const scalar searchRadiusSqr,
    List<pointIndexHit>& info
) const
{
    const PtrList<indexedOctree<treeDataEdge>>& edgeTrees = edgeTreesByType();

    DynamicList<pointIndexHit> dynEdgeHit(edgeTrees.size()*3);

    // Loop over all the feature edge types
    forAll(edgeTrees, i)
    {
        const auto& tree = edgeTrees[i];
        const auto& treeData = tree.shapes();

        // Pick up all the edges that intersect the search sphere
        labelList elems = tree.findSphere(sample, searchRadiusSqr);

        for (const label index : elems)
        {
            pointHit hitPoint = treeData.line(index).nearestDist(sample);

            // The index returned by indexedOctree is local its slice of
            // edges. Return the index into the complete edge list

            const label hitIndex = treeData.objectIndex(index);

            dynEdgeHit.append(pointIndexHit(hitPoint, hitIndex));
        }
    }

    info.transfer(dynEdgeHit);
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::extendedEdgeMesh::pointTree() const
{
    if (!pointTree_)
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4, ROOTVSMALL)
        );

        const labelList featurePointLabels = identity(nonFeatureStart_);

        pointTree_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(points(), featurePointLabels),

                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return *pointTree_;
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::extendedEdgeMesh::edgeTree() const
{
    if (!edgeTree_)
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4, ROOTVSMALL)
        );

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge(edges(), points()),  // All edges

                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return *edgeTree_;
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge>>&
Foam::extendedEdgeMesh::edgeTreesByType() const
{
    if (edgeTreesByType_.empty())
    {
        Random rndGen(872141);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4, ROOTVSMALL)
        );

        List<labelRange> sliceEdges(nEdgeTypes);

        // External edges
        sliceEdges[0].reset(externalStart_, (internalStart_ - externalStart_));

        // Internal edges
        sliceEdges[1].reset(internalStart_, (flatStart_ - internalStart_));

        // Flat edges
        sliceEdges[2].reset(flatStart_, (openStart_ - flatStart_));

        // Open edges
        sliceEdges[3].reset(openStart_, (multipleStart_ - openStart_));

        // Multiple edges
        sliceEdges[4].reset(multipleStart_, (edges().size() - multipleStart_));


        edgeTreesByType_.resize(nEdgeTypes);

        forAll(edgeTreesByType_, i)
        {
            edgeTreesByType_.set
            (
                i,
                new indexedOctree<treeDataEdge>
                (
                    // Selected edges
                    treeDataEdge(edges(), points(), sliceEdges[i]),

                    bb,     // bb
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }

    return edgeTreesByType_;
}


void Foam::extendedEdgeMesh::transfer(extendedEdgeMesh& mesh)
{
    if (&mesh == this)
    {
        return;  // Self-transfer is a no-op
    }

    edgeMesh::transfer(mesh);

    concaveStart_ = mesh.concaveStart_;
    mixedStart_ = mesh.mixedStart_;
    nonFeatureStart_ = mesh.nonFeatureStart_;
    internalStart_ = mesh.internalStart_;
    flatStart_ = mesh.flatStart_;
    openStart_ = mesh.openStart_;
    multipleStart_ = mesh.multipleStart_;
    normals_.transfer(mesh.normals_);
    normalVolumeTypes_.transfer(mesh.normalVolumeTypes_);
    edgeDirections_.transfer(mesh.edgeDirections_);
    normalDirections_.transfer(mesh.normalDirections_);
    edgeNormals_.transfer(mesh.edgeNormals_);
    featurePointNormals_.transfer(mesh.featurePointNormals_);
    featurePointEdges_.transfer(mesh.featurePointEdges_);
    regionEdges_.transfer(mesh.regionEdges_);
    pointTree_ = std::move(mesh.pointTree_);
    edgeTree_ = std::move(mesh.edgeTree_);
    edgeTreesByType_.transfer(mesh.edgeTreesByType_);

    mesh.clear();
}


void Foam::extendedEdgeMesh::clear()
{
    edgeMesh::clear();
    concaveStart_ = 0;
    mixedStart_ = 0;
    nonFeatureStart_ = 0;
    internalStart_ = 0;
    flatStart_ = 0;
    openStart_ = 0;
    multipleStart_ = 0;
    normals_.clear();
    normalVolumeTypes_.clear();
    edgeDirections_.clear();
    normalDirections_.clear();
    edgeNormals_.clear();
    featurePointNormals_.clear();
    featurePointEdges_.clear();
    regionEdges_.clear();
    pointTree_.reset(nullptr);
    edgeTree_.reset(nullptr);
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::add(const extendedEdgeMesh& fem)
{
    // Points
    // ~~~~~~

    // From current points into combined points
    labelList reversePointMap(points().size());
    labelList reverseFemPointMap(fem.points().size());

    label newPointi = 0;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = 0; i < fem.concaveStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Concave
    label newConcaveStart = newPointi;
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.concaveStart(); i < fem.mixedStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Mixed
    label newMixedStart = newPointi;
    for (label i = mixedStart(); i < nonFeatureStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.mixedStart(); i < fem.nonFeatureStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Non-feature
    label newNonFeatureStart = newPointi;
    for (label i = nonFeatureStart(); i < points().size(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.nonFeatureStart(); i < fem.points().size(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    pointField newPoints(newPointi);
    newPoints.rmap(points(), reversePointMap);
    newPoints.rmap(fem.points(), reverseFemPointMap);


    // Edges
    // ~~~~~

    // From current edges into combined edges
    labelList reverseEdgeMap(edges().size());
    labelList reverseFemEdgeMap(fem.edges().size());

    // External
    label newEdgeI = 0;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = 0; i < fem.internalStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Internal
    label newInternalStart = newEdgeI;
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.internalStart(); i < fem.flatStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Flat
    label newFlatStart = newEdgeI;
    for (label i = flatStart(); i < openStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.flatStart(); i < fem.openStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Open
    label newOpenStart = newEdgeI;
    for (label i = openStart(); i < multipleStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.openStart(); i < fem.multipleStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Multiple
    label newMultipleStart = newEdgeI;
    for (label i = multipleStart(); i < edges().size(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.multipleStart(); i < fem.edges().size(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    edgeList newEdges(newEdgeI);
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }
    forAll(fem.edges(), i)
    {
        const edge& e = fem.edges()[i];
        newEdges[reverseFemEdgeMap[i]] = edge
        (
            reverseFemPointMap[e[0]],
            reverseFemPointMap[e[1]]
        );
    }

    pointField newEdgeDirections
    (
        edgeDirections().size()
      + fem.edgeDirections().size()
    );
    newEdgeDirections.rmap(edgeDirections(), reverseEdgeMap);
    newEdgeDirections.rmap(fem.edgeDirections(), reverseFemEdgeMap);


    // Normals
    // ~~~~~~~

    // Combine normals
    DynamicField<point> newNormals
    (
        normals().size()
      + fem.normals().size()
    );
    newNormals.append(normals());
    newNormals.append(fem.normals());


    // Combine and re-index into newNormals
    labelListList newEdgeNormals
    (
        edgeNormals().size()
      + fem.edgeNormals().size()
    );

    UIndirectList<labelList>
    (
        newEdgeNormals,
        SubList<label>(reverseEdgeMap, edgeNormals().size())
    ) = edgeNormals();
    UIndirectList<labelList>
    (
        newEdgeNormals,
        SubList<label>(reverseFemEdgeMap, fem.edgeNormals().size())
    ) = fem.edgeNormals();

    forAll(fem.edgeNormals(), i)
    {
        const label mapI = reverseFemEdgeMap[i];
        labelList& en = newEdgeNormals[mapI];
        forAll(en, j)
        {
            en[j] += normals().size();
        }
    }


    // Combine and re-index into newFeaturePointNormals
    labelListList newFeaturePointNormals
    (
       featurePointNormals().size()
     + fem.featurePointNormals().size()
    );

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reverseFemPointMap, fem.featurePointNormals().size())
    ) = fem.featurePointNormals();

    forAll(fem.featurePointNormals(), i)
    {
        const label mapI = reverseFemPointMap[i];
        labelList& fn = newFeaturePointNormals[mapI];
        forAll(fn, j)
        {
            fn[j] += normals().size();
        }
    }


    // Combine regionEdges
    DynamicList<label> newRegionEdges
    (
        regionEdges().size()
      + fem.regionEdges().size()
    );
    forAll(regionEdges(), i)
    {
        newRegionEdges.append(reverseEdgeMap[regionEdges()[i]]);
    }
    forAll(fem.regionEdges(), i)
    {
        newRegionEdges.append(reverseFemEdgeMap[fem.regionEdges()[i]]);
    }


    // Assign
    // ~~~~~~

    // Transfer
    concaveStart_ = newConcaveStart;
    mixedStart_ = newMixedStart;
    nonFeatureStart_ = newNonFeatureStart;

    // Reset points and edges
    {
        edgeMesh newmesh(std::move(newPoints), std::move(newEdges));
        edgeMesh::transfer(newmesh);
    }

    // Transfer
    internalStart_ = newInternalStart;
    flatStart_ = newFlatStart;
    openStart_ = newOpenStart;
    multipleStart_ = newMultipleStart;

    edgeDirections_.transfer(newEdgeDirections);

    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);

    regionEdges_.transfer(newRegionEdges);

    pointTree_.reset(nullptr);
    edgeTree_.reset(nullptr);
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::flipNormals()
{
    // Points
    // ~~~~~~

    // From current points into new points
    labelList reversePointMap(identity(points().size()));

    // Flip convex and concave points

    label newPointi = 0;
    // Concave points become convex
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    // Convex points become concave
    label newConcaveStart = newPointi;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }


    // Edges
    // ~~~~~~

    // From current edges into new edges
    labelList reverseEdgeMap(identity(edges().size()));

    // Flip external and internal edges

    label newEdgeI = 0;
    // Internal become external
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    // External become internal
    label newInternalStart = newEdgeI;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }


    pointField newPoints(points().size());
    newPoints.rmap(points(), reversePointMap);

    edgeList newEdges(edges().size());
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }


    // Normals are flipped
    // ~~~~~~~~~~~~~~~~~~~

    pointField newEdgeDirections(edges().size());
    newEdgeDirections.rmap(-1.0*edgeDirections(), reverseEdgeMap);

    pointField newNormals(-1.0*normals());

    labelListList newEdgeNormals(edgeNormals().size());
    UIndirectList<labelList>(newEdgeNormals, reverseEdgeMap) = edgeNormals();

    labelListList newFeaturePointNormals(featurePointNormals().size());

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();

    labelList newRegionEdges(regionEdges().size());
    forAll(regionEdges(), i)
    {
        newRegionEdges[i] = reverseEdgeMap[regionEdges()[i]];
    }

    // Transfer
    concaveStart_ = newConcaveStart;

    // Reset points and edges
    {
        edgeMesh newmesh(std::move(newPoints), std::move(newEdges));
        edgeMesh::transfer(newmesh);
    }

    // Transfer
    internalStart_ = newInternalStart;

    edgeDirections_.transfer(newEdgeDirections);
    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);
    regionEdges_.transfer(newRegionEdges);

    pointTree_.reset(nullptr);
    edgeTree_.reset(nullptr);
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::autoMap
(
    const pointField& subPoints,
    const edgeList& subEdges,
    const labelList& pointMap,
    const labelList& edgeMap
)
{
    // Determine slicing for subEdges
    label subIntStart = edgeMap.size();
    label subFlatStart = edgeMap.size();
    label subOpenStart = edgeMap.size();
    label subMultipleStart = edgeMap.size();

    forAll(edgeMap, subEdgeI)
    {
        label edgeI = edgeMap[subEdgeI];
        if (edgeI >= internalStart() && subIntStart == edgeMap.size())
        {
            subIntStart = subEdgeI;
        }
        if (edgeI >= flatStart() && subFlatStart == edgeMap.size())
        {
            subFlatStart = subEdgeI;
        }
        if (edgeI >= openStart() && subOpenStart == edgeMap.size())
        {
            subOpenStart = subEdgeI;
        }
        if (edgeI >= multipleStart() && subMultipleStart == edgeMap.size())
        {
            subMultipleStart = subEdgeI;
        }
    }


    // Determine slicing for subPoints

    label subConcaveStart = pointMap.size();
    label subMixedStart = pointMap.size();
    label subNonFeatStart = pointMap.size();

    forAll(pointMap, subPointI)
    {
        label pointI = pointMap[subPointI];
        if (pointI >= concaveStart() && subConcaveStart == pointMap.size())
        {
            subConcaveStart = subPointI;
        }
        if (pointI >= mixedStart() && subMixedStart == pointMap.size())
        {
            subMixedStart = subPointI;
        }
        if
        (
            pointI >= nonFeatureStart()
         && subNonFeatStart == pointMap.size()
        )
        {
            subNonFeatStart = subPointI;
        }
    }



    // Compact region edges
    labelList subRegionEdges;
    {
        bitSet isRegionEdge(edges().size(), regionEdges());

        DynamicList<label> newRegionEdges(regionEdges().size());
        forAll(edgeMap, subEdgeI)
        {
            if (isRegionEdge.test(edgeMap[subEdgeI]))
            {
                newRegionEdges.append(subEdgeI);
            }
        }
        subRegionEdges.transfer(newRegionEdges);
    }


    labelListList subFeaturePointEdges;
    if (featurePointEdges().size())
    {
        subFeaturePointEdges.setSize(subNonFeatStart);
        for (label subPointI = 0; subPointI < subNonFeatStart; subPointI++)
        {
            label pointI = pointMap[subPointI];
            const labelList& pEdges = featurePointEdges()[pointI];

            labelList& subPEdges = subFeaturePointEdges[subPointI];
            subPEdges.setSize(pEdges.size());

            if (pEdges.size())
            {
                forAll(pEdges, i)
                {
                    subPEdges[i] = edgeMap[pEdges[i]];
                }
            }
        }
    }


    vectorField subEdgeDirections(edgeDirections(), edgeMap);

    // Find used normals
    labelList reverseNormalMap(normals().size(), -1);
    DynamicList<label> normalMap(normals().size());

    {
        bitSet isSubNormal(normals().size());
        for (label subPointI = 0; subPointI < subNonFeatStart; subPointI++)
        {
            label pointI = pointMap[subPointI];
            const labelList& pNormals = featurePointNormals()[pointI];

            isSubNormal.set(pNormals);
        }
        forAll(edgeMap, subEdgeI)
        {
            label edgeI = edgeMap[subEdgeI];
            const labelList& eNormals = edgeNormals()[edgeI];

            isSubNormal.set(eNormals);
        }

        forAll(isSubNormal, normalI)
        {
            if (isSubNormal.test(normalI))
            {
                label subNormalI = normalMap.size();
                reverseNormalMap[normalI] = subNormalI;
                normalMap.append(subNormalI);
            }
        }
    }


    // Use compaction map on data referencing normals
    labelListList subNormalDirections;

    if (normalDirections().size())
    {
        subNormalDirections.setSize(edgeMap.size());

        forAll(edgeMap, subEdgeI)
        {
            label edgeI = edgeMap[subEdgeI];
            const labelList& eNormals = normalDirections()[edgeI];

            labelList& subNormals = subNormalDirections[subEdgeI];
            subNormals.setSize(eNormals.size());
            forAll(eNormals, i)
            {
                if (eNormals[i] >= 0)
                {
                    subNormals[i] = reverseNormalMap[eNormals[i]];
                }
                else
                {
                    subNormals[i] = -1;
                }
            }
        }
    }

    labelListList subEdgeNormals(edgeMap.size());
    forAll(edgeMap, subEdgeI)
    {
        label edgeI = edgeMap[subEdgeI];
        const labelList& eNormals = edgeNormals()[edgeI];
        labelList& subNormals = subEdgeNormals[subEdgeI];

        subNormals = labelUIndList(reverseNormalMap, eNormals);
    }

    labelListList subPointNormals(pointMap.size());
    for (label subPointI = 0; subPointI < subNonFeatStart; subPointI++)
    {
        label pointI = pointMap[subPointI];
        const labelList& pNormals = featurePointNormals()[pointI];
        labelList& subNormals = subPointNormals[subPointI];

        subNormals = labelUIndList(reverseNormalMap, pNormals);
    }

    // Use compaction map to compact normal data
    vectorField subNormals(normals(), normalMap);

    List<extendedEdgeMesh::sideVolumeType> subNormalVolumeTypes;
    if (normalVolumeTypes().size())
    {
        subNormalVolumeTypes =
            UIndirectList<extendedEdgeMesh::sideVolumeType>
            (
                normalVolumeTypes(),
                normalMap
            );
    }

    extendedEdgeMesh subMesh
    (
        subPoints,
        subEdges,

        // Feature points slices
        subConcaveStart,
        subMixedStart,
        subNonFeatStart,

        // Feature edges slices
        subIntStart,
        subFlatStart,
        subOpenStart,
        subMultipleStart,

        // All normals
        subNormals,
        subNormalVolumeTypes,

        // Per edge edge vector
        subEdgeDirections,

        // Per edge list of normal indices
        subNormalDirections,
        // Per edge list of normal indices
        subEdgeNormals,

        // Per point list of normal indices
        subPointNormals,
        subFeaturePointEdges,
        subRegionEdges
    );

    transfer(subMesh);
}


void Foam::extendedEdgeMesh::trim
(
    const searchableSurface& surf,
    const volumeType volType,
    labelList& pointMap,
    labelList& edgeMap
)
{
    // Cut edges with the other surfaces

    labelList allPointMap;      // from all to original point
    labelList allEdgeMap;       // from all to original edge

    labelList pointsFromEdge;   // list of new points created by cutting
    labelList oldEdge;          // for each of these points the original edge
    labelList surfTri;          // for each of these points the surface triangle
    cut
    (
        surf,
        allPointMap,
        allEdgeMap,
        pointsFromEdge,
        oldEdge,
        surfTri
    );

    const label nOldPoints = points().size();

    // Remove outside edges and compact

    labelList subPointMap;  // sub to old points
    labelList subEdgeMap;   // sub to old edges
    select(surf, volType, subPointMap, subEdgeMap);

    // Update overall point maps
    pointMap = labelUIndList(allPointMap, subPointMap);
    edgeMap = labelUIndList(allEdgeMap, subEdgeMap);

    // Extract current point and edge status
    List<edgeStatus> edgeStat(edges().size());
    List<pointStatus> pointStat(points().size());
    forAll(edgeStat, edgeI)
    {
        edgeStat[edgeI] = getEdgeStatus(edgeI);
    }
    forAll(pointStat, pointI)
    {
        pointStat[pointI] = getPointStatus(pointI);
    }

    // Re-classify exposed points (from cutting)
    labelList oldPointToIndex(nOldPoints, -1);
    forAll(pointsFromEdge, i)
    {
        oldPointToIndex[pointsFromEdge[i]] = i;
    }
    forAll(subPointMap, pointI)
    {
        label oldPointI = subPointMap[pointI];
        label index = oldPointToIndex[oldPointI];
        if (index != -1)
        {
            pointStat[pointI] = classifyFeaturePoint(pointI);
        }
    }

    // Reset based on new point and edge status
    labelList sortedToOriginalPoint;
    labelList sortedToOriginalEdge;
    setFromStatus
    (
        pointStat,
        edgeStat,
        sortedToOriginalPoint,
        sortedToOriginalEdge
    );

    // Update the overall pointMap, edgeMap
    pointMap = labelUIndList(pointMap, sortedToOriginalPoint)();
    edgeMap = labelUIndList(edgeMap, sortedToOriginalEdge)();
}


void Foam::extendedEdgeMesh::setFromStatus
(
    const List<extendedEdgeMesh::pointStatus>& pointStat,
    const List<extendedEdgeMesh::edgeStatus>& edgeStat,
    labelList& sortedToOriginalPoint,
    labelList& sortedToOriginalEdge
)
{
    // Use pointStatus and edgeStatus to determine new ordering
    label pointConcaveStart;
    label pointMixedStart;
    label pointNonFeatStart;

    label edgeInternalStart;
    label edgeFlatStart;
    label edgeOpenStart;
    label edgeMultipleStart;
    sortedOrder
    (
        pointStat,
        edgeStat,
        sortedToOriginalPoint,
        sortedToOriginalEdge,

        pointConcaveStart,
        pointMixedStart,
        pointNonFeatStart,

        edgeInternalStart,
        edgeFlatStart,
        edgeOpenStart,
        edgeMultipleStart
    );


    // Order points and edges
    labelList reversePointMap(points().size(), -1);
    forAll(sortedToOriginalPoint, sortedI)
    {
        reversePointMap[sortedToOriginalPoint[sortedI]] = sortedI;
    }

    edgeList sortedEdges(UIndirectList<edge>(edges(), sortedToOriginalEdge)());
    forAll(sortedEdges, sortedI)
    {
        inplaceRenumber(reversePointMap, sortedEdges[sortedI]);
    }

    // Update local data
    autoMap
    (
        pointField(points(), sortedToOriginalPoint),
        sortedEdges,
        sortedToOriginalPoint,
        sortedToOriginalEdge
    );

    // Reset the slice starts
    concaveStart_ = pointConcaveStart;
    mixedStart_ = pointMixedStart;
    nonFeatureStart_ = pointNonFeatStart;
    internalStart_ = edgeInternalStart;
    flatStart_ = edgeFlatStart;
    openStart_ = edgeOpenStart;
    multipleStart_ = edgeMultipleStart;
}


bool Foam::extendedEdgeMesh::mergePointsAndSort
(
    const scalar mergeDist,
    labelList& pointMap,
    labelList& edgeMap
)
{
    const label nOldPoints = points().size();

    // Detect and merge collocated feature points
    labelList oldToMerged;
    label nNewPoints = Foam::mergePoints
    (
        points(),
        SMALL,
        false,
        oldToMerged
    );

    pointMap.setSize(nNewPoints);
    pointMap = -1;
    forAll(oldToMerged, oldI)
    {
        label newI = oldToMerged[oldI];
        if (pointMap[newI] == -1)
        {
            pointMap[newI] = oldI;
        }
    }

    // Renumber edges
    edgeList newEdges(edges().size());
    forAll(edges(), edgeI)
    {
        const edge& oldE = edges()[edgeI];
        newEdges[edgeI] = edge(oldToMerged[oldE[0]], oldToMerged[oldE[1]]);
    }

    // Shuffle basic information (reorders point data)
    autoMap
    (
        pointField(points(), pointMap),
        newEdges,
        pointMap,
        identity(newEdges.size())
    );

    // Re-classify the merged points
    List<edgeStatus> edgeStat(edges().size());
    forAll(edgeStat, edgeI)
    {
        edgeStat[edgeI] = getEdgeStatus(edgeI);
    }

    List<pointStatus> pointStat(points().size());
    forAll(pointStat, pointI)
    {
        pointStat[pointI] = getPointStatus(pointI);
    }

    // Re-classify merged points
    labelList nPoints(nNewPoints, Zero);
    forAll(oldToMerged, oldPointI)
    {
        nPoints[oldToMerged[oldPointI]]++;
    }

    forAll(nPoints, pointI)
    {
        if (nPoints[pointI] != 1)
        {
            pointStat[pointI] = classifyFeaturePoint(pointI);
        }
    }

    labelList sortedToOriginalPoint;
    setFromStatus
    (
        pointStat,
        edgeStat,
        sortedToOriginalPoint,
        edgeMap             // point merging above did not affect edge order
    );
    pointMap = labelUIndList(pointMap, sortedToOriginalPoint)();

    return nNewPoints != nOldPoints;
}


void Foam::extendedEdgeMesh::writeObj(const fileName& prefix) const
{
    Info<< nl << "Writing extendedEdgeMesh components to " << prefix
        << endl;

    edgeMesh::write(prefix + "_edgeMesh.obj");

    {
        OBJstream convexFtPtStr(prefix + "_convexFeaturePts.obj");
        Info<< "Writing " << concaveStart_
            << " convex feature points to " << convexFtPtStr.name() << endl;

        for (label i = 0; i < concaveStart_; i++)
        {
            convexFtPtStr.write(points()[i]);
        }
    }

    {
        OBJstream concaveFtPtStr(prefix + "_concaveFeaturePts.obj");
        Info<< "Writing " << mixedStart_-concaveStart_
            << " concave feature points to "
            << concaveFtPtStr.name() << endl;

        for (label i = concaveStart_; i < mixedStart_; i++)
        {
            concaveFtPtStr.write(points()[i]);
        }
    }

    {
        OBJstream mixedFtPtStr(prefix + "_mixedFeaturePts.obj");
        Info<< "Writing " << nonFeatureStart_-mixedStart_
            << " mixed feature points to " << mixedFtPtStr.name() << endl;

        for (label i = mixedStart_; i < nonFeatureStart_; i++)
        {
            mixedFtPtStr.write(points()[i]);
        }
    }

    {
        OBJstream mixedFtPtStructureStr(prefix+"_mixedFeaturePtsStructure.obj");
        Info<< "Writing "
            << nonFeatureStart_-mixedStart_
            << " mixed feature point structure to "
            << mixedFtPtStructureStr.name() << endl;

        for (label i = mixedStart_; i < nonFeatureStart_; i++)
        {
            const labelList& ptEds = pointEdges()[i];

            for (const label edgei : ptEds)
            {
                const edge& e = edges()[edgei];
                mixedFtPtStructureStr.write(e, points());
            }
        }
    }

    {
        OBJstream externalStr(prefix + "_externalEdges.obj");
        Info<< "Writing " << internalStart_-externalStart_
            << " external edges to " << externalStr.name() << endl;

        for (label i = externalStart_; i < internalStart_; i++)
        {
            const edge& e = edges()[i];
            externalStr.write(e, points());
        }
    }

    {
        OBJstream internalStr(prefix + "_internalEdges.obj");
        Info<< "Writing " << flatStart_-internalStart_
            << " internal edges to " << internalStr.name() << endl;

        for (label i = internalStart_; i < flatStart_; i++)
        {
            const edge& e = edges()[i];
            internalStr.write(e, points());
        }
    }

    {
        OBJstream flatStr(prefix + "_flatEdges.obj");
        Info<< "Writing " << openStart_-flatStart_
            << " flat edges to " << flatStr.name() << endl;

        for (label i = flatStart_; i < openStart_; i++)
        {
            const edge& e = edges()[i];
            flatStr.write(e, points());
        }
    }

    {
        OBJstream openStr(prefix + "_openEdges.obj");
        Info<< "Writing " << multipleStart_-openStart_
            << " open edges to " << openStr.name() << endl;

        for (label i = openStart_; i < multipleStart_; i++)
        {
            const edge& e = edges()[i];
            openStr.write(e, points());
        }
    }

    {
        OBJstream multipleStr(prefix + "_multipleEdges.obj");
        Info<< "Writing " << edges().size()-multipleStart_
            << " multiple edges to " << multipleStr.name() << endl;

        for (label i = multipleStart_; i < edges().size(); i++)
        {
            const edge& e = edges()[i];
            multipleStr.write(e, points());
        }
    }

    {
        OBJstream regionStr(prefix + "_regionEdges.obj");
        Info<< "Writing " << regionEdges_.size()
            << " region edges to " << regionStr.name() << endl;

        for (const label edgei : regionEdges_)
        {
            const edge& e = edges()[edgei];
            regionStr.write(e, points());
        }
    }

    {
        OBJstream edgeDirsStr(prefix + "_edgeDirections.obj");
        Info<< "Writing " << edgeDirections_.size()
            << " edge directions to " << edgeDirsStr.name() << endl;

        forAll(edgeDirections_, i)
        {
            const vector& eVec = edgeDirections_[i];
            const edge& e = edges()[i];

            edgeDirsStr.writeLine
            (
                points()[e.start()],
                points()[e.start()] + eVec
            );
        }
    }
}


void Foam::extendedEdgeMesh::writeStats(Ostream& os) const
{
    edgeMesh::writeStats(os);

    os  << indent << "point classification :" << nl;
    os  << incrIndent;
    os  << indent << "convex feature points          : "
        << setw(8) << concaveStart_-convexStart_
        //<< setw(8) << convexStart_
        << nl;
    os  << indent << "concave feature points         : "
        << setw(8) << mixedStart_-concaveStart_
        //<< setw(8) << concaveStart_
        << nl;
    os  << indent << "mixed feature points           : "
        << setw(8) << nonFeatureStart_-mixedStart_
        //<< setw(8) << mixedStart_
        << nl;
    os  << indent << "other (non-feature) points     : "
        << setw(8) << points().size()-nonFeatureStart_
        //<< setw(8) << nonFeatureStart_
        << nl;
    os  << decrIndent;

    os  << indent << "edge classification :" << nl;
    os  << incrIndent;
    os  << indent << "external (convex angle) edges  : "
        << setw(8) << internalStart_-externalStart_
        //<< setw(8) << externalStart_
        << nl;
    os  << indent << "internal (concave angle) edges : "
        << setw(8) << flatStart_-internalStart_
        //<< setw(8) << internalStart_
        << nl;
    os  << indent << "flat region edges              : "
        << setw(8) << openStart_-flatStart_
        //<< setw(8) << flatStart_
        << nl;
    os  << indent << "open edges                     : "
        << setw(8) << multipleStart_-openStart_
        //<< setw(8) << openStart_
        << nl;
    os  << indent << "multiply connected edges       : "
        << setw(8) << edges().size()-multipleStart_
        //<< setw(8) << multipleStart_
        << nl;
    os  << decrIndent;
}


Foam::extendedEdgeMesh::edgeStatus
Foam::extendedEdgeMesh::classifyEdge
(
    const List<vector>& norms,
    const labelList& edNorms,
    const vector& fC0tofC1
)
{
    label nEdNorms = edNorms.size();

    if (nEdNorms == 1)
    {
        return OPEN;
    }
    else if (nEdNorms == 2)
    {
        const vector& n0(norms[edNorms[0]]);
        const vector& n1(norms[edNorms[1]]);

        if ((n0 & n1) > cosNormalAngleTol_)
        {
            return FLAT;
        }
        else if ((fC0tofC1 & n0) > 0.0)
        {
            return INTERNAL;
        }
        else
        {
            return EXTERNAL;
        }
    }
    else if (nEdNorms > 2)
    {
        return MULTIPLE;
    }

    // There is a problem - the edge has no normals
    return NONE;
}


void Foam::extendedEdgeMesh::sortedOrder
(
    const List<extendedEdgeMesh::pointStatus>& pointStat,
    const List<extendedEdgeMesh::edgeStatus>& edgeStat,
    labelList& sortedToOriginalPoint,
    labelList& sortedToOriginalEdge,

    label& pointConcaveStart,
    label& pointMixedStart,
    label& pointNonFeatStart,

    label& edgeInternalStart,
    label& edgeFlatStart,
    label& edgeOpenStart,
    label& edgeMultipleStart
)
{
    sortedToOriginalPoint.setSize(pointStat.size());
    sortedToOriginalPoint = -1;

    sortedToOriginalEdge.setSize(edgeStat.size());
    sortedToOriginalEdge = -1;


    // Order edges
    // ~~~~~~~~~~~

    label nConvex = 0;
    label nConcave = 0;
    label nMixed = 0;
    // label nNonFeat = 0;

    forAll(pointStat, pointI)
    {
        switch (pointStat[pointI])
        {
            case extendedEdgeMesh::CONVEX:
                nConvex++;
            break;

            case extendedEdgeMesh::CONCAVE:
                nConcave++;
            break;

            case extendedEdgeMesh::MIXED:
                nMixed++;
            break;

            case extendedEdgeMesh::NONFEATURE:
                // ++nNonFeat;
            break;

            default:
                FatalErrorInFunction
                    << "Problem" << exit(FatalError);
            break;
        }
    }

    label convexStart = 0;
    label concaveStart = nConvex;
    label mixedStart = concaveStart+nConcave;
    label nonFeatStart = mixedStart+nMixed;


    // Copy to parameters
    pointConcaveStart = concaveStart;
    pointMixedStart = mixedStart;
    pointNonFeatStart = nonFeatStart;

    forAll(pointStat, pointI)
    {
        switch (pointStat[pointI])
        {
            case extendedEdgeMesh::CONVEX:
                sortedToOriginalPoint[convexStart++] = pointI;
            break;

            case extendedEdgeMesh::CONCAVE:
                sortedToOriginalPoint[concaveStart++] = pointI;
            break;

            case extendedEdgeMesh::MIXED:
                sortedToOriginalPoint[mixedStart++] = pointI;
            break;

            case extendedEdgeMesh::NONFEATURE:
                sortedToOriginalPoint[nonFeatStart++] = pointI;
            break;
        }
    }


    // Order edges
    // ~~~~~~~~~~~

    label nExternal = 0;
    label nInternal = 0;
    label nFlat = 0;
    label nOpen = 0;
    // label nMultiple = 0;

    forAll(edgeStat, edgeI)
    {
        switch (edgeStat[edgeI])
        {
            case extendedEdgeMesh::EXTERNAL:
                nExternal++;
            break;

            case extendedEdgeMesh::INTERNAL:
                nInternal++;
            break;

            case extendedEdgeMesh::FLAT:
                nFlat++;
            break;

            case extendedEdgeMesh::OPEN:
                nOpen++;
            break;

            case extendedEdgeMesh::MULTIPLE:
                // ++nMultiple;
            break;

            case extendedEdgeMesh::NONE:
            default:
                FatalErrorInFunction
                    << "Problem" << exit(FatalError);
            break;
        }
    }

    label externalStart = 0;
    label internalStart = nExternal;
    label flatStart = internalStart + nInternal;
    label openStart = flatStart + nFlat;
    label multipleStart = openStart + nOpen;


    // Copy to parameters
    edgeInternalStart = internalStart;
    edgeFlatStart = flatStart;
    edgeOpenStart = openStart;
    edgeMultipleStart = multipleStart;

    forAll(edgeStat, edgeI)
    {
        switch (edgeStat[edgeI])
        {
            case extendedEdgeMesh::EXTERNAL:
                sortedToOriginalEdge[externalStart++] = edgeI;
            break;

            case extendedEdgeMesh::INTERNAL:
                sortedToOriginalEdge[internalStart++] = edgeI;
            break;

            case extendedEdgeMesh::FLAT:
                sortedToOriginalEdge[flatStart++] = edgeI;
            break;

            case extendedEdgeMesh::OPEN:
                sortedToOriginalEdge[openStart++] = edgeI;
            break;

            case extendedEdgeMesh::MULTIPLE:
                sortedToOriginalEdge[multipleStart++] = edgeI;
            break;

            case extendedEdgeMesh::NONE:
            default:
                FatalErrorInFunction
                    << "Problem" << exit(FatalError);
            break;
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    Foam::extendedEdgeMesh::sideVolumeType& vt
)
{
    label type;
    is  >> type;

    vt = static_cast<Foam::extendedEdgeMesh::sideVolumeType>(type);

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Foam::extendedEdgeMesh::sideVolumeType& vt
)
{
    os  << static_cast<label>(vt);

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const extendedEdgeMesh& em)
{
    //fileFormats::extendedEdgeMeshFormat::write(os, em.points_, em.edges_);
    os  << "// points" << nl
        << em.points() << nl
        << "// edges" << nl
        << em.edges() << nl
        << "// concaveStart mixedStart nonFeatureStart" << nl
        << em.concaveStart_ << token::SPACE
        << em.mixedStart_ << token::SPACE
        << em.nonFeatureStart_ << nl
        << "// internalStart flatStart openStart multipleStart" << nl
        << em.internalStart_ << token::SPACE
        << em.flatStart_ << token::SPACE
        << em.openStart_ << token::SPACE
        << em.multipleStart_ << nl
        << "// normals" << nl
        << em.normals_ << nl
        << "// normal volume types" << nl
        << em.normalVolumeTypes_ << nl
        << "// normalDirections" << nl
        << em.normalDirections_ << nl
        << "// edgeNormals" << nl
        << em.edgeNormals_ << nl
        << "// featurePointNormals" << nl
        << em.featurePointNormals_ << nl
        << "// featurePointEdges" << nl
        << em.featurePointEdges_ << nl
        << "// regionEdges" << nl
        << em.regionEdges_
        << endl;

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, extendedEdgeMesh& em)
{
    //fileFormats::extendedEdgeMeshFormat::read(is, em.points_, em.edges_);
    is  >> static_cast<edgeMesh&>(em)
        >> em.concaveStart_
        >> em.mixedStart_
        >> em.nonFeatureStart_
        >> em.internalStart_
        >> em.flatStart_
        >> em.openStart_
        >> em.multipleStart_
        >> em.normals_
        >> em.normalVolumeTypes_
        >> em.normalDirections_
        >> em.edgeNormals_
        >> em.featurePointNormals_
        >> em.featurePointEdges_
        >> em.regionEdges_;

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
