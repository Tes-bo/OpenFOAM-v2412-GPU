/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

Application
    surfaceInflate

Group
    grpSurfaceUtilities

Description
    Inflates surface. WIP. Checks for overlaps and locally lowers
    inflation distance

Usage
    - surfaceInflate [OPTION]

    \param -checkSelfIntersection \n
    Includes checks for self-intersection.

    \param -nSmooth
    Specify number of smoothing iterations

    \param -featureAngle
    Specify a feature angle


    E.g. inflate surface by 20mm with 1.5 safety factor:
        surfaceInflate DTC-scaled.obj 0.02 1.5 -featureAngle 45 -nSmooth 2

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurfaceFields.H"
#include "triSurfaceMesh.H"
#include "triSurfaceGeoMesh.H"
#include "bitSet.H"
#include "OBJstream.H"
#include "surfaceFeatures.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calcVertexNormalWeight
(
    const triFace& f,
    const label pI,
    const vector& fN,
    const pointField& points
)
{
    label index = f.find(pI);

    if (index == -1)
    {
        FatalErrorInFunction
            << "Point not in face" << abort(FatalError);
    }

    const vector e1 = points[f[index]] - points[f[f.fcIndex(index)]];
    const vector e2 = points[f[index]] - points[f[f.rcIndex(index)]];

    return mag(fN)/(magSqr(e1)*magSqr(e2) + VSMALL);
}


tmp<vectorField> calcVertexNormals(const triSurface& surf)
{
    // Weighted average of normals of faces attached to the vertex
    // Weight = fA / (mag(e0)^2 * mag(e1)^2);
    auto tpointNormals = tmp<vectorField>::New(surf.nPoints(), Zero);
    auto& pointNormals = tpointNormals.ref();

    const pointField& points = surf.points();
    const labelListList& pointFaces = surf.pointFaces();
    const labelList& meshPoints = surf.meshPoints();

    forAll(pointFaces, pI)
    {
        const labelList& pFaces = pointFaces[pI];

        forAll(pFaces, fI)
        {
            const label faceI = pFaces[fI];
            const triFace& f = surf[faceI];

            vector areaNorm = f.areaNormal(points);

            scalar weight = calcVertexNormalWeight
            (
                f,
                meshPoints[pI],
                areaNorm,
                points
            );

            pointNormals[pI] += weight * areaNorm;
        }

        pointNormals[pI].normalise();
    }

    return tpointNormals;
}


tmp<vectorField> calcPointNormals
(
    const triSurface& s,
    const bitSet& isFeaturePoint,
    const List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    //const pointField pointNormals(s.pointNormals());
    tmp<vectorField> tpointNormals(calcVertexNormals(s));
    vectorField& pointNormals = tpointNormals.ref();


    // feature edges: create edge normals from edgeFaces only.
    {
        const labelListList& edgeFaces = s.edgeFaces();

        labelList nNormals(s.nPoints(), Zero);
        forAll(edgeStat, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const edge& e = s.edges()[edgeI];
                forAll(e, i)
                {
                    if (!isFeaturePoint.test(e[i]))
                    {
                        pointNormals[e[i]] = Zero;
                    }
                }
            }
        }

        forAll(edgeStat, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                // Get average edge normal
                vector n = Zero;
                for (const label facei : eFaces)
                {
                    n += s.faceNormals()[facei];
                }
                n /= eFaces.size();


                // Sum to points
                const edge& e = s.edges()[edgeI];
                forAll(e, i)
                {
                    if (!isFeaturePoint.test(e[i]))
                    {
                        pointNormals[e[i]] += n;
                        nNormals[e[i]]++;
                    }
                }
            }
        }

        forAll(nNormals, pointI)
        {
            if (nNormals[pointI] > 0)
            {
                pointNormals[pointI].normalise();
            }
        }
    }


    forAll(pointNormals, pointI)
    {
        if (mag(mag(pointNormals[pointI])-1) > SMALL)
        {
            FatalErrorInFunction
                << "unitialised"
                << exit(FatalError);
        }
    }

    return tpointNormals;
}


void detectSelfIntersections
(
    const triSurfaceMesh& s,
    bitSet& isEdgeIntersecting
)
{
    const edgeList& edges = s.edges();
    const indexedOctree<treeDataTriSurface>& tree = s.tree();
    const labelList& meshPoints = s.meshPoints();
    const tmp<pointField> tpoints(s.points());
    const pointField& points = tpoints();

    isEdgeIntersecting.resize_nocopy(edges.size());
    isEdgeIntersecting = false;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        pointIndexHit hitInfo
        (
            tree.findLine
            (
                points[meshPoints[e[0]]],
                points[meshPoints[e[1]]],
                treeDataTriSurface::findSelfIntersectOp(tree, edgeI)
            )
        );

        if (hitInfo.hit())
        {
            isEdgeIntersecting.set(edgeI);
        }
    }
}


label detectIntersectionPoints
(
    const scalar tolerance,
    const triSurfaceMesh& s,

    const scalar extendFactor,
    const pointField& initialPoints,
    const vectorField& displacement,

    const bool checkSelfIntersect,
    const bitSet& initialIsEdgeIntersecting,

    bitSet& isPointOnHitEdge,
    scalarField& scale
)
{
    const pointField initialLocalPoints(initialPoints, s.meshPoints());
    const List<labelledTri>& localFaces = s.localFaces();


    label nHits = 0;
    isPointOnHitEdge.setSize(s.nPoints());
    isPointOnHitEdge = false;


    // 1. Extrusion offset vectors intersecting new surface location
    {
        scalar tol = max(tolerance, 10*s.tolerance());

        // Collect all the edge vectors. Slightly shorten the edges to prevent
        // finding lots of intersections. The fast triangle intersection routine
        // has problems with rays passing through a point of the
        // triangle.
        pointField start(initialLocalPoints+tol*displacement);
        pointField end(initialLocalPoints+extendFactor*displacement);

        List<pointIndexHit> hits;
        s.findLineAny(start, end, hits);

        forAll(hits, pointI)
        {
            if
            (
                hits[pointI].hit()
            &&  !localFaces[hits[pointI].index()].found(pointI)
            )
            {
                scale[pointI] = max(0.0, scale[pointI]-0.2);

                isPointOnHitEdge.set(pointI);
                nHits++;
            }
        }
    }


    // 2. (new) surface self intersections
    if (checkSelfIntersect)
    {
        bitSet isEdgeIntersecting;
        detectSelfIntersections(s, isEdgeIntersecting);

        const edgeList& edges = s.edges();
        const tmp<pointField> tpoints(s.points());
        const pointField& points = tpoints();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            if (isEdgeIntersecting[edgeI] && !initialIsEdgeIntersecting[edgeI])
            {
                if (isPointOnHitEdge.set(e[0]))
                {
                    label start = s.meshPoints()[e[0]];
                    const point& pt = points[start];

                    Pout<< "Additional self intersection at "
                        << pt
                        << endl;

                    scale[e[0]] = max(0.0, scale[e[0]]-0.2);
                    nHits++;
                }
                if (isPointOnHitEdge.set(e[1]))
                {
                    label end = s.meshPoints()[e[1]];
                    const point& pt = points[end];

                    Pout<< "Additional self intersection at "
                        << pt
                        << endl;

                    scale[e[1]] = max(0.0, scale[e[1]]-0.2);
                    nHits++;
                }
            }
        }
    }

    return nHits;
}


tmp<scalarField> avg
(
    const triSurface& s,
    const scalarField& fld,
    const scalarField& edgeWeights
)
{
    auto tres = tmp<scalarField>::New(s.nPoints(), Zero);
    auto& res = tres.ref();

    scalarField sumWeight(s.nPoints(), Zero);

    const edgeList& edges = s.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const scalar w = edgeWeights[edgeI];

        res[e[0]] += w*fld[e[1]];
        sumWeight[e[0]] += w;

        res[e[1]] += w*fld[e[0]];
        sumWeight[e[1]] += w;
    }

    // Average
    // ~~~~~~~

    forAll(res, pointI)
    {
        if (mag(sumWeight[pointI]) < VSMALL)
        {
            // Unconnected point. Take over original value
            res[pointI] = fld[pointI];
        }
        else
        {
            res[pointI] /= sumWeight[pointI];
        }
    }

    return tres;
}


void minSmooth
(
    const triSurface& s,
    const bitSet& isAffectedPoint,
    const scalarField& fld,
    scalarField& newFld
)
{

    const edgeList& edges = s.edges();
    const pointField& points = s.points();
    const labelList& mp = s.meshPoints();

    scalarField edgeWeights(edges.size());
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        scalar w = mag(points[mp[e[0]]]-points[mp[e[1]]]);

        edgeWeights[edgeI] = 1.0/(max(w, SMALL));
    }

    tmp<scalarField> tavgFld = avg(s, fld, edgeWeights);

    const scalarField& avgFld = tavgFld();

    forAll(fld, pointI)
    {
        if (isAffectedPoint.test(pointI))
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
                //avgFld[pointI]
            );
        }
    }
}


void lloydsSmoothing
(
    const label nSmooth,
    triSurface& s,
    const bitSet& isFeaturePoint,
    const List<surfaceFeatures::edgeStatus>& edgeStat,
    const bitSet& isAffectedPoint
)
{
    const labelList& meshPoints = s.meshPoints();
    const edgeList& edges = s.edges();


    bitSet isSmoothPoint(isAffectedPoint);
    // Extend isSmoothPoint
    {
        bitSet newIsSmoothPoint(isSmoothPoint);
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];
            if (isSmoothPoint.test(e[0]))
            {
                newIsSmoothPoint.set(e[1]);
            }
            if (isSmoothPoint.test(e[1]))
            {
                newIsSmoothPoint.set(e[0]);
            }
        }
        Info<< "From points-to-smooth " << isSmoothPoint.count()
            << " to " << newIsSmoothPoint.count()
            << endl;
        isSmoothPoint.transfer(newIsSmoothPoint);
    }

    // Do some smoothing (Lloyds algorithm) around problematic points
    for (label i = 0; i < nSmooth; i++)
    {
        const labelListList& pointFaces = s.pointFaces();
        const pointField& faceCentres = s.faceCentres();

        pointField newPoints(s.points());
        forAll(isSmoothPoint, pointI)
        {
            if (isSmoothPoint[pointI] && !isFeaturePoint[pointI])
            {
                const labelList& pFaces = pointFaces[pointI];

                point avg(Zero);
                forAll(pFaces, pFaceI)
                {
                    avg += faceCentres[pFaces[pFaceI]];
                }
                avg /= pFaces.size();
                newPoints[meshPoints[pointI]] = avg;
            }
        }


        // Move points on feature edges only according to feature edges.

        const pointField& points = s.points();

        vectorField pointSum(s.nPoints(), Zero);
        labelList nPointSum(s.nPoints(), Zero);

        forAll(edges, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const edge& e = edges[edgeI];
                const point& start = points[meshPoints[e[0]]];
                const point& end = points[meshPoints[e[1]]];
                point eMid = 0.5*(start+end);
                pointSum[e[0]] += eMid;
                nPointSum[e[0]]++;
                pointSum[e[1]] += eMid;
                nPointSum[e[1]]++;
            }
        }

        forAll(pointSum, pointI)
        {
            if
            (
                isSmoothPoint[pointI]
             && isFeaturePoint[pointI]
             && nPointSum[pointI] >= 2
            )
            {
                newPoints[meshPoints[pointI]] =
                    pointSum[pointI]/nPointSum[pointI];
            }
        }


        s.movePoints(newPoints);


        // Extend isSmoothPoint
        {
            bitSet newIsSmoothPoint(isSmoothPoint);
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                if (isSmoothPoint[e[0]])
                {
                    newIsSmoothPoint.set(e[1]);
                }
                if (isSmoothPoint[e[1]])
                {
                    newIsSmoothPoint.set(e[0]);
                }
            }
            Info<< "From points-to-smooth " << isSmoothPoint.count()
                << " to " << newIsSmoothPoint.count()
                << endl;
            isSmoothPoint.transfer(newIsSmoothPoint);
        }
    }
}



// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Inflates surface according to point normals."
    );

    argList::noParallel();
    argList::addNote
    (
        "Creates inflated version of surface using point normals. "
        "Takes surface, distance to inflate and additional safety factor"
    );
    argList::addBoolOption
    (
        "checkSelfIntersection",
        "Also check for self-intersection"
    );
    argList::addOption
    (
        "nSmooth",
        "integer",
        "Number of smoothing iterations (default 20)"
    );
    argList::addOption
    (
        "featureAngle",
        "scalar",
        "Feature angle"
    );
    argList::addBoolOption
    (
        "debug",
        "Switch on additional debug information"
    );

    argList::addArgument("input", "The input surface file");
    argList::addArgument("distance", "The inflate distance");
    argList::addArgument("factor", "The extend safety factor [1,10]");

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    const auto inputName = args.get<word>(1);
    const auto distance = args.get<scalar>(2);
    const auto extendFactor = args.get<scalar>(3);
    const bool checkSelfIntersect = args.found("checkSelfIntersection");
    const auto nSmooth = args.getOrDefault<label>("nSmooth", 10);
    const auto featureAngle = args.getOrDefault<scalar>("featureAngle", 180);
    const bool debug = args.found("debug");


    Info<< "Inflating surface " << inputName
        << " according to point normals" << nl
        << "    distance               : " << distance << nl
        << "    safety factor          : " << extendFactor << nl
        << "    self intersection test : " << checkSelfIntersect << nl
        << "    smoothing iterations   : " << nSmooth << nl
        << "    feature angle          : " << featureAngle << nl
        << endl;


    if (extendFactor < 1 || extendFactor > 10)
    {
        FatalErrorInFunction
            << "Illegal safety factor " << extendFactor
            << ". It is usually 1..2"
            << exit(FatalError);
    }



    // Load triSurfaceMesh
    triSurfaceMesh s
    (
        IOobject
        (
            inputName,                          // name
            runTime.constant(),                 // instance
            "triSurface",                       // local
            runTime,                            // registry
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Mark features
    const surfaceFeatures features(s, 180.0-featureAngle);

    Info<< "Detected features:" << nl
        << "    feature points : " << features.featurePoints().size()
        << " out of " << s.nPoints() << nl
        << "    feature edges : " << features.featureEdges().size()
        << " out of " << s.nEdges() << nl
        << endl;

    bitSet isFeaturePoint(s.nPoints(), features.featurePoints());

    const List<surfaceFeatures::edgeStatus> edgeStat(features.toStatus());


    const pointField initialPoints(s.points());


    // Construct scale
    Info<< "Constructing field scale\n" << endl;
    triSurfacePointScalarField scale
    (
        IOobject
        (
            "scale",                            // name
            runTime.timeName(),                 // instance
            s,                                  // registry
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        s,
        scalar(1),
        dimLength
    );


    // Construct unit normals

    Info<< "Calculating vertex normals\n" << endl;
    const pointField pointNormals
    (
        calcPointNormals
        (
            s,
            isFeaturePoint,
            edgeStat
        )
    );


    // Construct pointDisplacement
    Info<< "Constructing field pointDisplacement\n" << endl;
    triSurfacePointVectorField pointDisplacement
    (
        IOobject
        (
            "pointDisplacement",                // name
            runTime.timeName(),                 // instance
            s,                                  // registry
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        s,
        dimLength,
        // - calculate with primitive fields (#2758)
        (distance*scale.field()*pointNormals)
    );


    const labelList& meshPoints = s.meshPoints();


    // Any point on any intersected edge in any of the iterations
    bitSet isScaledPoint(s.nPoints());


    // Detect any self intersections on initial mesh
    bitSet initialIsEdgeIntersecting;
    if (checkSelfIntersect)
    {
        detectSelfIntersections(s, initialIsEdgeIntersecting);
    }
    else
    {
        // Mark all edges as already self intersecting so avoid detecting any
        // new ones
        initialIsEdgeIntersecting.setSize(s.nEdges(), true);
    }


    // Inflate
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Move to new location
        pointField newPoints(initialPoints);
        forAll(meshPoints, i)
        {
            newPoints[meshPoints[i]] += pointDisplacement[i];
        }
        s.movePoints(newPoints);
        Info<< "Bounding box : " << s.bounds() << endl;


        // Work out scale from pointDisplacement
        forAll(scale, pointI)
        {
            if (s.pointFaces()[pointI].size())
            {
                scale[pointI] = mag(pointDisplacement[pointI])/distance;
            }
            else
            {
                scale[pointI] = 1.0;
            }
        }


        Info<< "Scale        : " << gAverage(scale) << endl;
        Info<< "Displacement : " << gAverage(pointDisplacement) << endl;



        // Detect any intersections and scale back
        bitSet isAffectedPoint;
        label nIntersections = detectIntersectionPoints
        (
            1e-9,       // intersection tolerance
            s,
            extendFactor,
            initialPoints,
            pointDisplacement,

            checkSelfIntersect,
            initialIsEdgeIntersecting,

            isAffectedPoint,
            scale
        );
        Info<< "Detected " << nIntersections << " intersections" << nl
            << endl;

        if (nIntersections == 0)
        {
            runTime.write();
            break;
        }


        // Accumulate all affected points
        forAll(isAffectedPoint, pointI)
        {
            if (isAffectedPoint.test(pointI))
            {
                isScaledPoint.set(pointI);
            }
        }

        // Smear out lowering of scale so any edges not found are
        // still included
        for (label i = 0; i < nSmooth; i++)
        {
            triSurfacePointScalarField oldScale(scale);
            oldScale.rename("oldScale");
            minSmooth
            (
                s,
                bitSet(s.nPoints(), true),
                oldScale,
                scale
            );
        }


        // From scale update the pointDisplacement
        // - calculate with primitive fields (#2758)
        pointDisplacement.normalise();
        pointDisplacement.field() *= distance*scale.field();


        // Do some smoothing (Lloyds algorithm)
        lloydsSmoothing(nSmooth, s, isFeaturePoint, edgeStat, isAffectedPoint);

        // Update pointDisplacement
        const tmp<pointField> tpoints(s.points());
        const pointField& pts = tpoints();

        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];
            pointDisplacement[i] = pts[meshPointI]-initialPoints[meshPointI];
        }


        // Write
        runTime.write();

        runTime.printExecutionTime(Info);
    }


    Info<< "Detected overall intersecting points " << isScaledPoint.count()
        << " out of " << s.nPoints() << nl << endl;


    if (debug)
    {
        OBJstream str(runTime.path()/"isScaledPoint.obj");

        for (const label pointi : isScaledPoint)
        {
            str.write(initialPoints[meshPoints[pointi]]);
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
