/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "surfaceFeatures.H"
#include "triSurface.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "treeDataPoint.H"
#include "meshTools.H"
#include "Fstream.H"
#include "unitConversion.H"
#include "edgeHashes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceFeatures, 0);

    const scalar surfaceFeatures::parallelTolerance = sin(degToRad(1.0));


//! \cond fileScope
//  Check if the point is on the line
static bool onLine(const Foam::point& p, const linePointRef& line)
{
    const point& a = line.start();
    const point& b = line.end();

    if
    (
        (p.x() < min(a.x(), b.x()) || p.x() > max(a.x(), b.x()))
     || (p.y() < min(a.y(), b.y()) || p.y() > max(a.y(), b.y()))
     || (p.z() < min(a.z(), b.z()) || p.z() > max(a.z(), b.z()))
    )
    {
        return false;
    }

    return true;
}
//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::surfaceFeatures::edgeNearest
(
    const linePointRef& line,
    const point& sample
)
{
    pointHit eHit = line.nearestDist(sample);

    // Classification of position on edge.
    label endPoint;

    if (eHit.hit())
    {
        // Nearest point is on edge itself.
        // Note: might be at or very close to endpoint. We should use tolerance
        // here.
        endPoint = -1;
    }
    else
    {
        // Nearest point has to be one of the end points. Find out
        // which one.
        if
        (
            eHit.point().distSqr(line.start())
          < eHit.point().distSqr(line.end())
        )
        {
            endPoint = 0;
        }
        else
        {
            endPoint = 1;
        }
    }

    return pointIndexHit(eHit, endPoint);
}


// Go from selected edges only to a value for every edge
Foam::List<Foam::surfaceFeatures::edgeStatus> Foam::surfaceFeatures::toStatus()
 const
{
    List<edgeStatus> edgeStat(surf_.nEdges(), NONE);

    // Region edges first
    for (label i = 0; i < externalStart_; i++)
    {
        edgeStat[featureEdges_[i]] = REGION;
    }

    // External edges
    for (label i = externalStart_; i < internalStart_; i++)
    {
        edgeStat[featureEdges_[i]] = EXTERNAL;
    }

    // Internal edges
    for (label i = internalStart_; i < featureEdges_.size(); i++)
    {
        edgeStat[featureEdges_[i]] = INTERNAL;
    }

    return edgeStat;
}


// Set from value for every edge
void Foam::surfaceFeatures::setFromStatus
(
    const List<edgeStatus>& edgeStat,
    const scalar includedAngle
)
{
    // Count

    label nRegion = 0;
    label nExternal = 0;
    label nInternal = 0;

    forAll(edgeStat, edgeI)
    {
        if (edgeStat[edgeI] == REGION)
        {
            nRegion++;
        }
        else if (edgeStat[edgeI] == EXTERNAL)
        {
            nExternal++;
        }
        else if (edgeStat[edgeI] == INTERNAL)
        {
            nInternal++;
        }
    }

    externalStart_ = nRegion;
    internalStart_ = externalStart_ + nExternal;


    // Copy

    featureEdges_.setSize(internalStart_ + nInternal);

    label regionI = 0;
    label externalI = externalStart_;
    label internalI = internalStart_;

    forAll(edgeStat, edgeI)
    {
        if (edgeStat[edgeI] == REGION)
        {
            featureEdges_[regionI++] = edgeI;
        }
        else if (edgeStat[edgeI] == EXTERNAL)
        {
            featureEdges_[externalI++] = edgeI;
        }
        else if (edgeStat[edgeI] == INTERNAL)
        {
            featureEdges_[internalI++] = edgeI;
        }
    }

    const scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    calcFeatPoints(edgeStat, minCos);
}


//construct feature points where more than 2 feature edges meet
void Foam::surfaceFeatures::calcFeatPoints
(
    const List<edgeStatus>& edgeStat,
    const scalar minCos
)
{
    DynamicList<label> featurePoints(surf_.nPoints()/1000);

    const labelListList& pointEdges = surf_.pointEdges();
    const edgeList& edges = surf_.edges();
    const pointField& localPoints = surf_.localPoints();

    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];

        label nFeatEdges = 0;

        forAll(pEdges, i)
        {
            if (edgeStat[pEdges[i]] != NONE)
            {
                nFeatEdges++;
            }
        }

        if (nFeatEdges > 2)
        {
            featurePoints.append(pointi);
        }
        else if (nFeatEdges == 2)
        {
            // Check the angle between the two edges
            DynamicList<vector> edgeVecs(2);

            forAll(pEdges, i)
            {
                const label edgeI = pEdges[i];

                if (edgeStat[edgeI] != NONE)
                {
                    vector vec = edges[edgeI].vec(localPoints);
                    scalar magVec = mag(vec);
                    if (magVec > SMALL)
                    {
                        edgeVecs.append(vec/magVec);
                    }
                }
            }

            if (edgeVecs.size() == 2 && mag(edgeVecs[0] & edgeVecs[1]) < minCos)
            {
                featurePoints.append(pointi);
            }
        }
    }

    featurePoints_.transfer(featurePoints);
}


void Foam::surfaceFeatures::classifyFeatureAngles
(
    const labelListList& edgeFaces,
    List<edgeStatus>& edgeStat,
    const scalar minCos,
    const bool geometricTestOnly
) const
{
    const vectorField& faceNormals = surf_.faceNormals();
    const pointField& points = surf_.points();

    // Special case: minCos=1
    bool selectAll = (mag(minCos-1.0) < SMALL);

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() != 2)
        {
            // Non-manifold. What to do here? Is region edge? external edge?
            edgeStat[edgeI] = REGION;
        }
        else
        {
            label face0 = eFaces[0];
            label face1 = eFaces[1];

            if
            (
                !geometricTestOnly
             && surf_[face0].region() != surf_[face1].region()
            )
            {
                edgeStat[edgeI] = REGION;
            }
            else if
            (
                selectAll
             || ((faceNormals[face0] & faceNormals[face1]) < minCos)
            )
            {
                // Check if convex or concave by looking at angle
                // between face centres and normal
                vector f0Tof1 =
                    surf_[face1].centre(points)
                  - surf_[face0].centre(points);

                if ((f0Tof1 & faceNormals[face0]) >= 0.0)
                {
                    edgeStat[edgeI] = INTERNAL;
                }
                else
                {
                    edgeStat[edgeI] = EXTERNAL;
                }
            }
        }
    }
}


// Returns next feature edge connected to pointi with correct value.
Foam::label Foam::surfaceFeatures::nextFeatEdge
(
    const List<edgeStatus>& edgeStat,
    const labelList& featVisited,
    const label unsetVal,
    const label prevEdgeI,
    const label vertI
) const
{
    const labelList& pEdges = surf_.pointEdges()[vertI];

    label nextEdgeI = -1;

    forAll(pEdges, i)
    {
        label edgeI = pEdges[i];

        if
        (
            edgeI != prevEdgeI
         && edgeStat[edgeI] != NONE
         && featVisited[edgeI] == unsetVal
        )
        {
            if (nextEdgeI == -1)
            {
                nextEdgeI = edgeI;
            }
            else
            {
                // More than one feature edge to choose from. End of segment.
                return -1;
            }
        }
    }

    return nextEdgeI;
}


// Finds connected feature edges by walking from prevEdgeI in direction of
// prevPointi. Marks feature edges visited in featVisited by assigning them
// the current feature line number. Returns cumulative length of edges walked.
// Works in one of two modes:
// - mark : step to edges with featVisited = -1.
//          Mark edges visited with currentFeatI.
// - clear : step to edges with featVisited = currentFeatI
//           Mark edges visited with -2 and erase from feature edges.
Foam::surfaceFeatures::labelScalar Foam::surfaceFeatures::walkSegment
(
    const bool mark,
    const List<edgeStatus>& edgeStat,
    const label startEdgeI,
    const label startPointi,
    const label currentFeatI,
    labelList& featVisited
)
{
    label edgeI = startEdgeI;

    label vertI = startPointi;

    scalar visitedLength = 0.0;

    label nVisited = 0;

    if (featurePoints_.found(startPointi))
    {
        // Do not walk across feature points

        return labelScalar(nVisited, visitedLength);
    }

    //
    // Now we have:
    //    edgeI : first edge on this segment
    //    vertI : one of the endpoints of this segment
    //
    // Start walking, marking off edges (in featVisited)
    // as we go along.
    //

    label unsetVal;

    if (mark)
    {
        unsetVal = -1;
    }
    else
    {
        unsetVal = currentFeatI;
    }

    do
    {
        // Step to next feature edge with value unsetVal
        edgeI = nextFeatEdge(edgeStat, featVisited, unsetVal, edgeI, vertI);

        if (edgeI == -1 || edgeI == startEdgeI)
        {
            break;
        }

        // Mark with current value. If in clean mode also remove feature edge.

        if (mark)
        {
            featVisited[edgeI] = currentFeatI;
        }
        else
        {
            featVisited[edgeI] = -2;
        }


        // Step to next vertex on edge
        const edge& e = surf_.edges()[edgeI];

        vertI = e.otherVertex(vertI);


        // Update cumulative length
        visitedLength += e.mag(surf_.localPoints());

        nVisited++;

        if (nVisited > surf_.nEdges())
        {
            Warning<< "walkSegment : reached iteration limit in walking "
                << "feature edges on surface from edge:" << startEdgeI
                << " vertex:" << startPointi << nl
                << "Returning with large length" << endl;

            return labelScalar(nVisited, GREAT);
        }
    }
    while (true);

    return labelScalar(nVisited, visitedLength);
}


//- Divide into multiple normal bins
//  - return REGION if != 2 normals
//  - return REGION if 2 normals that make feature angle
//  - otherwise return NONE and set normals,bins
//
// This has been relocated from surfaceFeatureExtract and could be cleaned
// up some more.
//
Foam::surfaceFeatures::edgeStatus
Foam::surfaceFeatures::surfaceFeatures::checkFlatRegionEdge
(
    const scalar tol,
    const scalar includedAngle,
    const label edgeI,
    const point& leftPoint
) const
{
    const triSurface& surf = surf_;

    const edge& e = surf.edges()[edgeI];
    const labelList& eFaces = surf.edgeFaces()[edgeI];

    // Bin according to normal

    DynamicList<vector> normals(2);
    DynamicList<labelList> bins(2);

    forAll(eFaces, eFacei)
    {
        const vector& n = surf.faceNormals()[eFaces[eFacei]];

        // Find the normal in normals
        label index = -1;
        forAll(normals, normalI)
        {
            if (mag(n & normals[normalI]) > (1-tol))
            {
                index = normalI;
                break;
            }
        }

        if (index != -1)
        {
            bins[index].append(eFacei);
        }
        else if (normals.size() >= 2)
        {
            // Would be third normal. Mark as feature.
            //Pout<< "** at edge:" << surf.localPoints()[e[0]]
            //    << surf.localPoints()[e[1]]
            //    << " have normals:" << normals
            //    << " and " << n << endl;
            return surfaceFeatures::REGION;
        }
        else
        {
            normals.append(n);
            bins.append(labelList(1, eFacei));
        }
    }


    // Check resulting number of bins
    if (bins.size() == 1)
    {
        // Note: should check here whether they are two sets of faces
        // that are planar or indeed 4 faces al coming together at an edge.
        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]]
        //    << " have single normal:" << normals[0]
        //    << endl;
        return surfaceFeatures::NONE;
    }
    else
    {
        // Two bins. Check if normals make an angle

        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]] << nl
        //    << "    normals:" << normals << nl
        //    << "    bins   :" << bins << nl
        //    << endl;

        if (includedAngle >= 0)
        {
            scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

            forAll(eFaces, i)
            {
                const vector& ni = surf.faceNormals()[eFaces[i]];
                for (label j=i+1; j<eFaces.size(); j++)
                {
                    const vector& nj = surf.faceNormals()[eFaces[j]];
                    if (mag(ni & nj) < minCos)
                    {
                        //Pout<< "have sharp feature between normal:" << ni
                        //    << " and " << nj << endl;

                        // Is feature. Keep as region or convert to
                        // feature angle? For now keep as region.
                        return surfaceFeatures::REGION;
                    }
                }
            }
        }

        // So now we have two normals bins but need to make sure both
        // bins have the same regions in it.

        // 1. store + or - region number depending
        //    on orientation of triangle in bins[0]
        const labelList& bin0 = bins[0];
        labelList regionAndNormal(bin0.size());
        forAll(bin0, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin0[i]]];
            const auto dir = t.edgeDirection(e);

            if (dir > 0)
            {
                regionAndNormal[i] = t.region()+1;
            }
            else if (dir == 0)
            {
                FatalErrorInFunction
                    << exit(FatalError);
            }
            else
            {
                regionAndNormal[i] = -(t.region()+1);
            }
        }

        // 2. check against bin1
        const labelList& bin1 = bins[1];
        labelList regionAndNormal1(bin1.size());
        forAll(bin1, i)
        {
            const labelledTri& t = surf.localFaces()[eFaces[bin1[i]]];
            const auto dir = t.edgeDirection(e);

            label myRegionAndNormal;
            if (dir > 0)
            {
                myRegionAndNormal = t.region()+1;
            }
            else
            {
                myRegionAndNormal = -(t.region()+1);
            }

            regionAndNormal1[i] = myRegionAndNormal;

            label index = regionAndNormal.find(-myRegionAndNormal);
            if (index == -1)
            {
                // Not found.
                //Pout<< "cannot find region " << myRegionAndNormal
                //    << " in regions " << regionAndNormal << endl;

                return surfaceFeatures::REGION;
            }
        }

        // Passed all checks, two normal bins with the same contents.
        //Pout<< "regionAndNormal:" << regionAndNormal << endl;
        //Pout<< "myRegionAndNormal:" << regionAndNormal1 << endl;
    }

    return surfaceFeatures::NONE;
}

/*

// TBD. Big problem for duplicate triangles with opposing normals: we
// don't know which one of the duplicates gets found on either side so
// the normal might be + or -. Hence on e.g. motorBike.obj we see feature
// lines where there shouldn't be
Foam::surfaceFeatures::edgeStatus
Foam::surfaceFeatures::surfaceFeatures::checkBinRegion
(
    const label edgei,
    const labelList& bin0,
    const labelList& bin1
) const
{
    const triSurface& surf = surf_;
    const labelList& eFaces = surf.edgeFaces()[edgei];
    const edge& e = surf.edges()[edgei];

    // 1. store + or - region number depending
    //    on orientation of triangle in bins[0]
    labelList regionAndNormal(bin0.size());
    forAll(bin0, i)
    {
        const labelledTri& t = surf.localFaces()[eFaces[bin0[i]]];
        const auto dir = t.edgeDirection(e);

        if (dir > 0)
        {
            regionAndNormal[i] = t.region()+1;
        }
        else if (dir == 0)
        {
            FatalErrorInFunction
                << exit(FatalError);
        }
        else
        {
            regionAndNormal[i] = -(t.region()+1);
        }
    }

    // 2. check against bin1
    labelList regionAndNormal1(bin1.size());
    forAll(bin1, i)
    {
        const labelledTri& t = surf.localFaces()[eFaces[bin1[i]]];
        const auto dir = t.edgeDirection(e);

        label myRegionAndNormal;
        if (dir > 0)
        {
            myRegionAndNormal = t.region()+1;
        }
        else if (dir == 0)
        {
            myRegionAndNormal = 0;
            FatalErrorInFunction
                << exit(FatalError);
        }
        else
        {
            myRegionAndNormal = -(t.region()+1);
        }

        regionAndNormal1[i] = myRegionAndNormal;

        label index = regionAndNormal.find(-myRegionAndNormal);
        if (index == -1)
        {
            // Not found.
            //Pout<< "cannot find region " << myRegionAndNormal
            //    << " in regions " << regionAndNormal << endl;

            return surfaceFeatures::REGION;
        }
    }
    return surfaceFeatures::NONE;
}


Foam::surfaceFeatures::edgeStatus
Foam::surfaceFeatures::surfaceFeatures::checkFlatRegionEdge
(
    const scalar tol,
    const scalar includedAngle,
    const label edgei,
    const point& leftPoint
) const
{
    const triSurface& surf = surf_;
    const edge& e = surf.edges()[edgei];
    const auto& mp = surf.meshPoints();
    const point eMid(edge(mp[e[0]], mp[e[1]]).centre(surf.points()));
    const labelList& eFaces = surf.edgeFaces()[edgei];

    vector leftVec(leftPoint-eMid);
    leftVec.normalise();

    // Bin according to normal and location w.r.t. first face
    plane pl(eMid, leftVec);


    DynamicList<labelList> leftBins(2);
    DynamicList<labelList> rightBins(2);
    DynamicList<vector> leftNormals(2);
    DynamicList<vector> rightNormals(2);

    // Append first face since this is what leftPoint was created from in the
    // first place
    leftNormals.append(surf.faceNormals()[eFaces[0]]);
    leftBins.append(labelList(1, 0));

    for (label eFacei = 1; eFacei < eFaces.size(); ++eFacei)
    {
        const label facei = eFaces[eFacei];

        const bool isLeft(pl.signedDistance(surf.faceCentres()[facei]) > 0);

        DynamicList<labelList>& bins = (isLeft ? leftBins : rightBins);
        DynamicList<vector>& normals = (isLeft ? leftNormals : rightNormals);

        const vector& n = surf.faceNormals()[facei];

        // Find the normal in normals
        label index = -1;
        forAll(normals, normalI)
        {
            if (mag(n & normals[normalI]) > (1-tol))
            {
                index = normalI;
                break;
            }
        }

        if (index != -1)
        {
            bins[index].append(eFacei);
//            Pout<< "edge:" << edgei << " verts:" << e
//                << " found existing normal bin:" << index
//                << " after:" << flatOutput(bins[index])
//                << endl;
        }
        else if ((leftNormals.size()+rightNormals.size()) >= 2)
        {
            // Would be third normal. Mark as feature.
            //Pout<< "** at edge:" << surf.localPoints()[e[0]]
            //    << surf.localPoints()[e[1]]
            //    << " have normals:" << normals
            //    << " and " << n << endl;
            return surfaceFeatures::REGION;
        }
        else
        {
            normals.append(n);
            bins.append(labelList(1, eFacei));
        }
    }

    // Check resulting number of bins
    if ((leftBins.size()+rightBins.size()) == 1)
    {
        // Note: should check here whether they are two sets of faces
        // that are planar or indeed 4 faces al coming together at an edge.
        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]]
        //    << " have single normal:" << normals[0]
        //    << endl;
        return surfaceFeatures::NONE;
    }
    else
    {
        // Two bins. Check if normals make an angle

        //Pout<< "** at edge:"
        //    << surf.localPoints()[e[0]]
        //    << surf.localPoints()[e[1]] << nl
        //    << "    normals:" << normals << nl
        //    << "    bins   :" << bins << nl
        //    << endl;

        if (includedAngle >= 0)
        {
            scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

            forAll(eFaces, i)
            {
                const vector& ni = surf.faceNormals()[eFaces[i]];
                for (label j=i+1; j<eFaces.size(); j++)
                {
                    const vector& nj = surf.faceNormals()[eFaces[j]];
                    if (mag(ni & nj) < minCos)
                    {
                        //Pout<< "have sharp feature between normal:" << ni
                        //    << " and " << nj << endl;

                        // Is feature. Keep as region or convert to
                        // feature angle? For now keep as region.
                        return surfaceFeatures::REGION;
                    }
                }
            }
        }

        // So now we have two normals bins but need to make sure both
        // bins have the same regions in it.

        if (leftBins.size() == 2)
        {
            return checkBinRegion
            (
                edgei,
                leftBins[0],
                leftBins[1]
            );
        }
        else if (leftBins.size() == 1 && rightBins.size() == 1)
        {
            return checkBinRegion
            (
                edgei,
                leftBins[0],
                rightBins[0]
            );
        }
        else if (rightBins.size() == 2)
        {
            return checkBinRegion
            (
                edgei,
                rightBins[0],
                rightBins[1]
            );
        }
        else
        {
            FatalErrorInFunction << "leftBins:" << leftBins
                << " rightBins:" << rightBins << exit(FatalError);
        }
    }

    return surfaceFeatures::NONE;
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeatures::surfaceFeatures(const triSurface& surf)
:
    surf_(surf),
    featurePoints_(0),
    featureEdges_(0),
    externalStart_(0),
    internalStart_(0)
{}


// Construct from components.
Foam::surfaceFeatures::surfaceFeatures
(
    const triSurface& surf,
    const labelList& featurePoints,
    const labelList& featureEdges,
    const label externalStart,
    const label internalStart
)
:
    surf_(surf),
    featurePoints_(featurePoints),
    featureEdges_(featureEdges),
    externalStart_(externalStart),
    internalStart_(externalStart)
{}


// Construct from surface, angle and min length
Foam::surfaceFeatures::surfaceFeatures
(
    const triSurface& surf,
    const scalar includedAngle,
    const scalar minLen,
    const label minElems,
    const bool geometricTestOnly
)
:
    surf_(surf),
    featurePoints_(0),
    featureEdges_(0),
    externalStart_(0),
    internalStart_(0)
{
    findFeatures(includedAngle, geometricTestOnly);

    if (minLen > 0 || minElems > 0)
    {
        trimFeatures(minLen, minElems, includedAngle);
    }
}


Foam::surfaceFeatures::surfaceFeatures
(
    const triSurface& surf,
    const dictionary& featInfoDict
)
:
    surf_(surf),
    featurePoints_(featInfoDict.lookup("featurePoints")),
    featureEdges_(featInfoDict.lookup("featureEdges")),
    externalStart_(featInfoDict.get<label>("externalStart")),
    internalStart_(featInfoDict.get<label>("internalStart"))
{}


Foam::surfaceFeatures::surfaceFeatures
(
    const triSurface& surf,
    const fileName& fName
)
:
    surf_(surf),
    featurePoints_(0),
    featureEdges_(0),
    externalStart_(0),
    internalStart_(0)
{
    IFstream str(fName);

    dictionary featInfoDict(str);

    featInfoDict.readEntry("featureEdges", featureEdges_);
    featInfoDict.readEntry("featurePoints", featurePoints_);
    featInfoDict.readEntry("externalStart", externalStart_);
    featInfoDict.readEntry("internalStart", internalStart_);
}


Foam::surfaceFeatures::surfaceFeatures
(
    const triSurface& surf,
    const pointField& points,
    const edgeList& edges,
    const scalar mergeTol,
    const bool geometricTestOnly
)
:
    surf_(surf),
    featurePoints_(0),
    featureEdges_(0),
    externalStart_(0),
    internalStart_(0)
{
    // Match edge mesh edges with the triSurface edges

    const labelListList& surfEdgeFaces = surf_.edgeFaces();
    const edgeList& surfEdges = surf_.edges();

    scalar mergeTolSqr = sqr(mergeTol);

    EdgeMap<label> dynFeatEdges(2*edges.size());
    DynamicList<labelList> dynFeatureEdgeFaces(edges.size());

    labelList edgeLabel;

    nearestFeatEdge
    (
        edges,
        points,
        mergeTolSqr,
        edgeLabel     // label of surface edge or -1
    );

    label count = 0;
    forAll(edgeLabel, sEdgeI)
    {
        const label sEdge = edgeLabel[sEdgeI];

        if (sEdge == -1)
        {
            continue;
        }

        dynFeatEdges.insert(surfEdges[sEdge], count++);
        dynFeatureEdgeFaces.append(surfEdgeFaces[sEdge]);
    }

    // Find whether an edge is external or internal
    List<edgeStatus> edgeStat(dynFeatEdges.size(), NONE);

    classifyFeatureAngles
    (
        dynFeatureEdgeFaces,
        edgeStat,
        GREAT,
        geometricTestOnly
    );

    // Transfer the edge status to a list encompassing all edges in the surface
    // so that calcFeatPoints can be used.
    List<edgeStatus> allEdgeStat(surf_.nEdges(), NONE);

    forAll(allEdgeStat, eI)
    {
        const auto iter = dynFeatEdges.cfind(surfEdges[eI]);

        if (iter.good())
        {
            allEdgeStat[eI] = edgeStat[iter.val()];
        }
    }

    edgeStat.clear();
    dynFeatEdges.clear();

    setFromStatus(allEdgeStat, GREAT);
}


Foam::surfaceFeatures::surfaceFeatures(const surfaceFeatures& sf)
:
    surf_(sf.surface()),
    featurePoints_(sf.featurePoints()),
    featureEdges_(sf.featureEdges()),
    externalStart_(sf.externalStart()),
    internalStart_(sf.internalStart())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::surfaceFeatures::selectFeatureEdges
(
    const bool regionEdges,
    const bool externalEdges,
    const bool internalEdges
) const
{
    DynamicList<label> selectedEdges;

    if (regionEdges)
    {
        selectedEdges.setCapacity(selectedEdges.size() + nRegionEdges());

        for (label i = 0; i < externalStart_; i++)
        {
            selectedEdges.append(featureEdges_[i]);
        }
    }

    if (externalEdges)
    {
        selectedEdges.setCapacity(selectedEdges.size() + nExternalEdges());

        for (label i = externalStart_; i < internalStart_; i++)
        {
            selectedEdges.append(featureEdges_[i]);
        }
    }

    if (internalEdges)
    {
        selectedEdges.setCapacity(selectedEdges.size() + nInternalEdges());

        for (label i = internalStart_; i < featureEdges_.size(); i++)
        {
            selectedEdges.append(featureEdges_[i]);
        }
    }

    return selectedEdges.shrink();
}


void Foam::surfaceFeatures::findFeatures
(
    const scalar includedAngle,
    const bool geometricTestOnly
)
{
    scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    // Per edge whether is feature edge.
    List<edgeStatus> edgeStat(surf_.nEdges(), NONE);

    classifyFeatureAngles
    (
        surf_.edgeFaces(),
        edgeStat,
        minCos,
        geometricTestOnly
    );

    setFromStatus(edgeStat, includedAngle);
}


// Remove small strings of edges. First assigns string number to
// every edge and then works out the length of them.
Foam::labelList Foam::surfaceFeatures::trimFeatures
(
    const scalar minLen,
    const label minElems,
    const scalar includedAngle
)
{
    // Convert feature edge list to status per edge.
    List<edgeStatus> edgeStat(toStatus());


    // Mark feature edges according to connected lines.
    // -1: unassigned
    // -2: part of too small a feature line
    // >0: feature line number
    labelList featLines(surf_.nEdges(), -1);

    // Current featureline number.
    label featI = 0;

    // Current starting edge
    label startEdgeI = 0;

    do
    {
        // Find unset featureline
        for (; startEdgeI < edgeStat.size(); startEdgeI++)
        {
            if
            (
                edgeStat[startEdgeI] != NONE
             && featLines[startEdgeI] == -1
            )
            {
                // Found unassigned feature edge.
                break;
            }
        }

        if (startEdgeI == edgeStat.size())
        {
            // No unset feature edge found. All feature edges have line number
            // assigned.
            break;
        }

        featLines[startEdgeI] = featI;

        const edge& startEdge = surf_.edges()[startEdgeI];

        // Walk 'left' and 'right' from startEdge.
        labelScalar leftPath =
            walkSegment
            (
                true,           // 'mark' mode
                edgeStat,
                startEdgeI,     // edge, not yet assigned to a featureLine
                startEdge[0],   // start of edge
                featI,          // Mark value
                featLines       // Mark for all edges
            );

        labelScalar rightPath =
            walkSegment
            (
                true,
                edgeStat,
                startEdgeI,
                startEdge[1],   // end of edge
                featI,
                featLines
            );

        if
        (
            (
                leftPath.len_
              + rightPath.len_
              + startEdge.mag(surf_.localPoints())
              < minLen
            )
         || (leftPath.n_ + rightPath.n_ + 1 < minElems)
        )
        {
            // Rewalk same route (recognizable by featLines == featI)
            // to reset featLines.

            featLines[startEdgeI] = -2;

            walkSegment
            (
                false,          // 'clean' mode
                edgeStat,
                startEdgeI,     // edge, not yet assigned to a featureLine
                startEdge[0],   // start of edge
                featI,          // Unset value
                featLines       // Mark for all edges
            );

            walkSegment
            (
                false,
                edgeStat,
                startEdgeI,
                startEdge[1],   // end of edge
                featI,
                featLines
            );
        }
        else
        {
            featI++;
        }
    }
    while (true);

    // Unmark all feature lines that have featLines=-2
    forAll(featureEdges_, i)
    {
        label edgeI = featureEdges_[i];

        if (featLines[edgeI] == -2)
        {
            edgeStat[edgeI] = NONE;
        }
    }

    // Convert back to edge labels
    setFromStatus(edgeStat, includedAngle);

    return featLines;
}


void Foam::surfaceFeatures::excludeBox
(
    List<edgeStatus>& edgeStat,
    const treeBoundBox& bb
) const
{
    deleteBox(edgeStat, bb, true);
}


void Foam::surfaceFeatures::subsetBox
(
    List<edgeStatus>& edgeStat,
    const treeBoundBox& bb
) const
{
    deleteBox(edgeStat, bb, false);
}


void Foam::surfaceFeatures::deleteBox
(
    List<edgeStatus>& edgeStat,
    const treeBoundBox& bb,
    const bool removeInside
) const
{
    const edgeList& surfEdges = surf_.edges();
    const pointField& surfLocalPoints = surf_.localPoints();

    forAll(edgeStat, edgei)
    {
        const point eMid = surfEdges[edgei].centre(surfLocalPoints);

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}


void Foam::surfaceFeatures::subsetPlane
(
    List<edgeStatus>& edgeStat,
    const plane& cutPlane
) const
{
    const edgeList& surfEdges = surf_.edges();
    const pointField& pts = surf_.points();
    const labelList& meshPoints = surf_.meshPoints();

    forAll(edgeStat, edgei)
    {
        const edge& e = surfEdges[edgei];

        const point& p0 = pts[meshPoints[e.start()]];
        const point& p1 = pts[meshPoints[e.end()]];
        const linePointRef line(p0, p1);

        // If edge does not intersect the plane, delete.
        scalar intersect = cutPlane.lineIntersect(line);

        point featPoint = intersect * (p1 - p0) + p0;

        if (!onLine(featPoint, line))
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}


void Foam::surfaceFeatures::excludeOpen
(
    List<edgeStatus>& edgeStat
) const
{
    forAll(edgeStat, edgei)
    {
        if (surf_.edgeFaces()[edgei].size() == 1)
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}


void Foam::surfaceFeatures::excludeNonManifold
(
    List<edgeStatus>& edgeStat
) const
{
    forAll(edgeStat, edgei)
    {
        if (surf_.edgeFaces()[edgei].size() > 2)
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}


//- Divide into multiple normal bins
//  - return REGION if != 2 normals
//  - return REGION if 2 normals that make feature angle
//  - otherwise return NONE and set normals,bins
void Foam::surfaceFeatures::checkFlatRegionEdge
(
    List<edgeStatus>& edgeStat,
    const scalar tol,
    const scalar includedAngle
) const
{
    forAll(edgeStat, edgei)
    {
        if (edgeStat[edgei] == surfaceFeatures::REGION)
        {
            const labelList& eFaces = surf_.edgeFaces()[edgei];

            if (eFaces.size() > 2 && (eFaces.size() % 2) == 0)
            {
                const point& leftPoint = surf_.faceCentres()[eFaces[0]];

                edgeStat[edgei] = checkFlatRegionEdge
                (
                    tol,
                    includedAngle,
                    edgei,
                    leftPoint
                );
            }
        }
    }
}


void Foam::surfaceFeatures::writeDict(Ostream& os) const
{
    dictionary featInfoDict;
    featInfoDict.add("externalStart", externalStart_);
    featInfoDict.add("internalStart", internalStart_);
    featInfoDict.add("featureEdges", featureEdges_);
    featInfoDict.add("featurePoints", featurePoints_);

    featInfoDict.write(os);
}


void Foam::surfaceFeatures::write(const fileName& fName) const
{
    OFstream os(fName);
    writeDict(os);
}


void Foam::surfaceFeatures::writeObj(const fileName& prefix) const
{
    OFstream regionStr(prefix + "_regionEdges.obj");
    Pout<< "Writing region edges to " << regionStr.name() << endl;

    label verti = 0;
    for (label i = 0; i < externalStart_; i++)
    {
        const edge& e = surf_.edges()[featureEdges_[i]];

        meshTools::writeOBJ(regionStr, surf_.localPoints()[e[0]]); verti++;
        meshTools::writeOBJ(regionStr, surf_.localPoints()[e[1]]); verti++;
        regionStr << "l " << verti-1 << ' ' << verti << endl;
    }


    OFstream externalStr(prefix + "_externalEdges.obj");
    Pout<< "Writing external edges to " << externalStr.name() << endl;

    verti = 0;
    for (label i = externalStart_; i < internalStart_; i++)
    {
        const edge& e = surf_.edges()[featureEdges_[i]];

        meshTools::writeOBJ(externalStr, surf_.localPoints()[e[0]]); verti++;
        meshTools::writeOBJ(externalStr, surf_.localPoints()[e[1]]); verti++;
        externalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream internalStr(prefix + "_internalEdges.obj");
    Pout<< "Writing internal edges to " << internalStr.name() << endl;

    verti = 0;
    for (label i = internalStart_; i < featureEdges_.size(); i++)
    {
        const edge& e = surf_.edges()[featureEdges_[i]];

        meshTools::writeOBJ(internalStr, surf_.localPoints()[e[0]]); verti++;
        meshTools::writeOBJ(internalStr, surf_.localPoints()[e[1]]); verti++;
        internalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream pointStr(prefix + "_points.obj");
    Pout<< "Writing feature points to " << pointStr.name() << endl;

    for (const label pointi : featurePoints_)
    {
        meshTools::writeOBJ(pointStr, surf_.localPoints()[pointi]);
    }
}


void Foam::surfaceFeatures::writeStats(Ostream& os) const
{
    os  << "Feature set:" << nl
        << "    points : " << this->featurePoints().size() << nl
        << "    edges  : " << this->featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << this->nRegionEdges() << nl
        << "        external edges : " << this->nExternalEdges() << nl
        << "        internal edges : " << this->nInternalEdges() << endl;
}


// Get nearest vertex on patch for every point of surf in pointSet.
Foam::Map<Foam::label> Foam::surfaceFeatures::nearestSamples
(
    const labelList& pointLabels,
    const pointField& samples,
    const scalarField& maxDistSqr
) const
{
    // Build tree out of all samples.

    // Define bound box here (gcc-4.8.5)
    const treeBoundBox overallBb(samples);

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(samples),
        overallBb,
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );
    const auto& treeData = ppTree.shapes();

    // From patch point to surface point
    Map<label> nearest(2*pointLabels.size());

    const pointField& surfPoints = surf_.localPoints();

    forAll(pointLabels, i)
    {
        const label surfPointi = pointLabels[i];

        const point& surfPt = surfPoints[surfPointi];

        pointIndexHit info = ppTree.findNearest
        (
            surfPt,
            maxDistSqr[i]
        );

        if (!info.hit())
        {
            FatalErrorInFunction
                << "Problem for point "
                << surfPointi << " in tree " << ppTree.bb()
                << abort(FatalError);
        }

        label sampleI = info.index();

        if (treeData.centre(sampleI).distSqr(surfPt) < maxDistSqr[sampleI])
        {
            nearest.insert(sampleI, surfPointi);
        }
    }


    if (debug)
    {
        //
        // Dump to obj file
        //

        Pout<< "Dumping nearest surface feature points to nearestSamples.obj"
            << endl;

        OFstream objStream("nearestSamples.obj");

        label vertI = 0;
        forAllConstIters(nearest, iter)
        {
            meshTools::writeOBJ(objStream, samples[iter.key()]); vertI++;
            meshTools::writeOBJ(objStream, surfPoints[iter.val()]); vertI++;
            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest sample point for regularly sampled points along
// selected edges. Return map from sample to edge label
Foam::Map<Foam::label> Foam::surfaceFeatures::nearestSamples
(
    const labelList& selectedEdges,
    const pointField& samples,
    const scalarField& sampleDist,
    const scalarField& maxDistSqr,
    const scalar minSampleDist
) const
{
    const pointField& surfPoints = surf_.localPoints();
    const edgeList& surfEdges = surf_.edges();

    scalar maxSearchSqr = max(maxDistSqr);

    // Define bound box here (gcc-4.8.5)
    const treeBoundBox overallBb(samples);

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(samples),
        overallBb,
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    // From patch point to surface edge.
    Map<label> nearest(2*selectedEdges.size());

    forAll(selectedEdges, i)
    {
        label surfEdgeI = selectedEdges[i];

        const edge& e = surfEdges[surfEdgeI];

        if (debug && (i % 1000) == 0)
        {
            Pout<< "looking at surface feature edge " << surfEdgeI
                << " verts:" << e << " points:" << surfPoints[e[0]]
                << ' ' << surfPoints[e[1]] << endl;
        }

        // Normalized edge vector
        vector eVec = e.vec(surfPoints);
        scalar eMag = mag(eVec);
        eVec /= eMag;


        //
        // Sample along edge
        //

        bool exit = false;

        // Coordinate along edge (0 .. eMag)
        scalar s = 0.0;

        while (true)
        {
            point edgePoint(surfPoints[e.start()] + s*eVec);

            pointIndexHit info = ppTree.findNearest
            (
                edgePoint,
                maxSearchSqr
            );

            if (!info.hit())
            {
                // No point close enough to surface edge.
                break;
            }

            label sampleI = info.index();

            if (info.point().distSqr(edgePoint) < maxDistSqr[sampleI])
            {
                nearest.insert(sampleI, surfEdgeI);
            }

            if (exit)
            {
                break;
            }

            // Step to next sample point using local distance.
            // Truncate to max 1/minSampleDist samples per feature edge.
            s += max(minSampleDist*eMag, sampleDist[sampleI]);

            if (s >= (1-minSampleDist)*eMag)
            {
                // Do one more sample, at edge endpoint
                s = eMag;
                exit = true;
            }
        }
    }



    if (debug)
    {
        // Dump to obj file

        Pout<< "Dumping nearest surface edges to nearestEdges.obj"
            << endl;

        OFstream objStream("nearestEdges.obj");

        label vertI = 0;
        forAllConstIters(nearest, iter)
        {
            const label sampleI = iter.key();

            const edge& e = surfEdges[iter.val()];

            meshTools::writeOBJ(objStream, samples[sampleI]); vertI++;

            point nearPt =
                e.line(surfPoints).nearestDist(samples[sampleI]).point();

            meshTools::writeOBJ(objStream, nearPt); vertI++;

            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest edge on patch for regularly sampled points along the feature
// edges. Return map from patch edge to feature edge.
//
// Q: using point-based sampleDist and maxDist (distance to look around
// each point). Should they be edge-based e.g. by averaging or max()?
Foam::Map<Foam::pointIndexHit> Foam::surfaceFeatures::nearestEdges
(
    const labelList& selectedEdges,
    const edgeList& sampleEdges,
    const labelList& selectedSampleEdges,
    const pointField& samplePoints,
    const scalarField& sampleDist,
    const scalarField& maxDistSqr,
    const scalar minSampleDist
) const
{
    // Build tree out of selected sample edges.
    indexedOctree<treeDataEdge> ppTree
    (
        treeDataEdge
        (
            sampleEdges,
            samplePoints,
            selectedSampleEdges
        ),
        treeBoundBox(samplePoints), // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    const pointField& surfPoints = surf_.localPoints();
    const edgeList& surfEdges = surf_.edges();

    scalar maxSearchSqr = max(maxDistSqr);

    Map<pointIndexHit> nearest(2*sampleEdges.size());

    //
    // Loop over all selected edges. Sample at regular intervals. Find nearest
    // sampleEdges (using octree)
    //

    forAll(selectedEdges, i)
    {
        label surfEdgeI = selectedEdges[i];

        const edge& e = surfEdges[surfEdgeI];

        if (debug && (i % 1000) == 0)
        {
            Pout<< "looking at surface feature edge " << surfEdgeI
                << " verts:" << e << " points:" << surfPoints[e[0]]
                << ' ' << surfPoints[e[1]] << endl;
        }

        // Normalized edge vector
        vector eVec = e.vec(surfPoints);
        scalar eMag = mag(eVec);
        eVec /= eMag;


        //
        // Sample along edge
        //

        bool exit = false;

        // Coordinate along edge (0 .. eMag)
        scalar s = 0.0;

        while (true)
        {
            point edgePoint(surfPoints[e.start()] + s*eVec);

            pointIndexHit info = ppTree.findNearest
            (
                edgePoint,
                maxSearchSqr
            );

            if (!info.hit())
            {
                // No edge close enough to surface edge.
                break;
            }

            label index = info.index();

            label sampleEdgeI = ppTree.shapes().objectIndex(index);

            const edge& e = sampleEdges[sampleEdgeI];

            if (info.point().distSqr(edgePoint) < maxDistSqr[e.start()])
            {
                nearest.insert
                (
                    sampleEdgeI,
                    pointIndexHit(true, edgePoint, surfEdgeI)
                );
            }

            if (exit)
            {
                break;
            }

            // Step to next sample point using local distance.
            // Truncate to max 1/minSampleDist samples per feature edge.
            // s += max(minSampleDist*eMag, sampleDist[e.start()]);
            s += 0.01*eMag;

            if (s >= (1-minSampleDist)*eMag)
            {
                // Do one more sample, at edge endpoint
                s = eMag;
                exit = true;
            }
        }
    }


    if (debug)
    {
        // Dump to obj file

        Pout<< "Dumping nearest surface feature edges to nearestEdges.obj"
            << endl;

        OFstream objStream("nearestEdges.obj");

        label vertI = 0;
        forAllConstIters(nearest, iter)
        {
            const label sampleEdgeI = iter.key();

            const edge& sampleEdge = sampleEdges[sampleEdgeI];

            // Write line from edgeMid to point on feature edge
            meshTools::writeOBJ(objStream, sampleEdge.centre(samplePoints));
            vertI++;

            meshTools::writeOBJ(objStream, iter.val().point());
            vertI++;

            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest surface edge for every sample. Return in form of
// labelLists giving surfaceEdge label&intersectionpoint.
void Foam::surfaceFeatures::nearestSurfEdge
(
    const labelList& selectedEdges,
    const pointField& samples,
    scalar searchSpanSqr,   // Search span
    labelList& edgeLabel,
    labelList& edgeEndPoint,
    pointField& edgePoint
) const
{
    edgeLabel.setSize(samples.size());
    edgeEndPoint.setSize(samples.size());
    edgePoint.setSize(samples.size());

    const pointField& localPoints = surf_.localPoints();

    treeBoundBox searchDomain(localPoints);
    searchDomain.inflate(0.1);

    indexedOctree<treeDataEdge> ppTree
    (
        treeDataEdge
        (
            surf_.edges(),
            localPoints,
            selectedEdges
        ),
        searchDomain,  // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );
    const auto& treeData = ppTree.shapes();

    forAll(samples, i)
    {
        const point& sample = samples[i];

        pointIndexHit info = ppTree.findNearest
        (
            sample,
            searchSpanSqr
        );

        if (!info.hit())
        {
            edgeLabel[i] = -1;
        }
        else
        {
            // Need to recalculate to classify edgeEndPoint
            pointIndexHit pHit = edgeNearest
            (
                treeData.line(info.index()),
                sample
            );

            edgeLabel[i] = treeData.objectIndex(info.index());
            edgePoint[i] = pHit.point();
            edgeEndPoint[i] = pHit.index();
        }
    }
}


// Get nearest point on nearest feature edge for every sample (is edge)
void Foam::surfaceFeatures::nearestSurfEdge
(
    const labelList& selectedEdges,

    const edgeList& sampleEdges,
    const labelList& selectedSampleEdges,
    const pointField& samplePoints,

    const vector& searchSpan,   // Search span
    labelList& edgeLabel,       // label of surface edge or -1
    pointField& pointOnEdge,    // point on above edge
    pointField& pointOnFeature  // point on sample edge
) const
{
    edgeLabel.setSize(selectedSampleEdges.size());
    pointOnEdge.setSize(selectedSampleEdges.size());
    pointOnFeature.setSize(selectedSampleEdges.size());

    treeBoundBox searchDomain(surf_.localPoints());

    indexedOctree<treeDataEdge> ppTree
    (
        treeDataEdge
        (
            surf_.edges(),
            surf_.localPoints(),
            selectedEdges
        ),
        searchDomain,   // overall search domain
        8,              // maxLevel
        10,             // leafsize
        3.0             // duplicity
    );
    const auto& treeData = ppTree.shapes();

    forAll(selectedSampleEdges, i)
    {
        const edge& e = sampleEdges[selectedSampleEdges[i]];

        linePointRef edgeLine = e.line(samplePoints);

        point eMid(edgeLine.centre());

        treeBoundBox tightest(eMid - searchSpan, eMid + searchSpan);

        pointIndexHit info = ppTree.findNearest
        (
            edgeLine,
            tightest,
            pointOnEdge[i]
        );

        if (!info.hit())
        {
            edgeLabel[i] = -1;
        }
        else
        {
            edgeLabel[i] = treeData.objectIndex(info.index());
            pointOnFeature[i] = info.point();
        }
    }
}


void Foam::surfaceFeatures::nearestFeatEdge
(
    const edgeList& edges,
    const pointField& points,
    scalar searchSpanSqr,   // Search span
    labelList& edgeLabel
) const
{
    edgeLabel = labelList(surf_.nEdges(), -1);

    treeBoundBox searchDomain(points);
    searchDomain.inflate(0.1);

    indexedOctree<treeDataEdge> ppTree
    (
        treeDataEdge(edges, points),  // All edges

        searchDomain,  // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );
    const auto& treeData = ppTree.shapes();

    const edgeList& surfEdges = surf_.edges();
    const pointField& surfLocalPoints = surf_.localPoints();

    forAll(surfEdges, edgeI)
    {
        const edge& sample = surfEdges[edgeI];

        const point& startPoint = surfLocalPoints[sample.start()];
        const point& midPoint = sample.centre(surfLocalPoints);

        pointIndexHit infoMid = ppTree.findNearest
        (
            midPoint,
            searchSpanSqr
        );

        if (infoMid.hit())
        {
            const vector surfEdgeDir = midPoint - startPoint;

            const vector featEdgeDir = treeData.line(infoMid.index()).vec();

            // Check that the edges are nearly parallel
            if (mag(surfEdgeDir ^ featEdgeDir) < parallelTolerance)
            {
                edgeLabel[edgeI] = edgeI;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::surfaceFeatures::operator=(const surfaceFeatures& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    if (&surf_ != &rhs.surface())
    {
        FatalErrorInFunction
            << "Operating on different surfaces"
            << abort(FatalError);
    }

    featurePoints_ = rhs.featurePoints();
    featureEdges_ = rhs.featureEdges();
    externalStart_ = rhs.externalStart();
    internalStart_ = rhs.internalStart();
}


// ************************************************************************* //
