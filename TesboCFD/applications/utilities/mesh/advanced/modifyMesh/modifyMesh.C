/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    modifyMesh

Group
    grpMeshAdvancedUtilities

Description
    Manipulate mesh elements.

    Actions are:
        (boundary)points:
            - move

        (boundary)edges:
            - split and move introduced point

        (boundary)faces:
            - split(triangulate) and move introduced point

        edges:
            - collapse

        cells:
            - split into polygonal base pyramids around newly introduced mid
              point

    Is a bit of a loose collection of mesh change drivers.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "boundaryCutter.H"
#include "cellSplitter.H"
#include "edgeCollapser.H"
#include "meshTools.H"
#include "Pair.H"
#include "globalIndex.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Locate point on patch. Returns (mesh) point label.
label findPoint(const primitivePatch& pp, const point& nearPoint)
{
    const pointField& points = pp.points();
    const labelList& meshPoints = pp.meshPoints();

    // Find nearest and next nearest
    scalar minDistSqr = GREAT;
    label minI = -1;

    scalar almostMinDistSqr = GREAT;
    label almostMinI = -1;

    for (const label pointi : meshPoints)
    {
        scalar distSqr = nearPoint.distSqr(points[pointi]);

        if (distSqr < minDistSqr)
        {
            almostMinDistSqr = minDistSqr;
            almostMinI = minI;

            minDistSqr = distSqr;
            minI = pointi;
        }
        else if (distSqr < almostMinDistSqr)
        {
            almostMinDistSqr = distSqr;
            almostMinI = pointi;
        }
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest point      : " << minI
        << " distance " <<  Foam::sqrt(minDistSqr)
        << " at " << points[minI] << nl
        << "    next nearest point : " << almostMinI
        << " distance " <<  Foam::sqrt(almostMinDistSqr)
        << " at " << points[almostMinI] << endl;

    if (almostMinDistSqr < 4*minDistSqr)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}


// Locate edge on patch. Return mesh edge label.
label findEdge
(
    const primitiveMesh& mesh,
    const primitivePatch& pp,
    const point& nearPoint
)
{
    const pointField& localPoints = pp.localPoints();
    const pointField& points = pp.points();
    const labelList& meshPoints = pp.meshPoints();
    const edgeList& edges = pp.edges();

    // Find nearest and next nearest
    scalar minDist = GREAT;
    label minI = -1;

    scalar almostMinDist = GREAT;
    label almostMinI = -1;

    for (const edge& e : edges)
    {
        pointHit pHit(e.line(localPoints).nearestDist(nearPoint));

        if (pHit.hit())
        {
            if (pHit.distance() < minDist)
            {
                almostMinDist = minDist;
                almostMinI = minI;

                minDist = pHit.distance();
                minI = meshTools::findEdge
                (
                    mesh,
                    meshPoints[e[0]],
                    meshPoints[e[1]]
                );
            }
            else if (pHit.distance() < almostMinDist)
            {
                almostMinDist = pHit.distance();
                almostMinI = meshTools::findEdge
                (
                    mesh,
                    meshPoints[e[0]],
                    meshPoints[e[1]]
                );
            }
        }
    }

    if (minI == -1)
    {
        Info<< "Did not find edge close to point " << nearPoint << endl;

        return -1;
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest edge      : " << minI
        << " distance " << minDist << " endpoints "
        << mesh.edges()[minI].line(points) << nl
        << "    next nearest edge : " << almostMinI
        << " distance " << almostMinDist << " endpoints "
        << mesh.edges()[almostMinI].line(points) << nl
        << endl;

    if (almostMinDist < 2*minDist)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}


// Find face on patch. Return mesh face label.
label findFace
(
    const primitiveMesh& mesh,
    const primitivePatch& pp,
    const point& nearPoint
)
{
    const pointField& points = pp.points();

    // Find nearest and next nearest
    scalar minDist = GREAT;
    label minI = -1;

    scalar almostMinDist = GREAT;
    label almostMinI = -1;

    forAll(pp, patchFacei)
    {
        pointHit pHit(pp[patchFacei].nearestPoint(nearPoint, points));

        if (pHit.hit())
        {
            if (pHit.distance() < minDist)
            {
                almostMinDist = minDist;
                almostMinI = minI;

                minDist = pHit.distance();
                minI = patchFacei + mesh.nInternalFaces();
            }
            else if (pHit.distance() < almostMinDist)
            {
                almostMinDist = pHit.distance();
                almostMinI = patchFacei + mesh.nInternalFaces();
            }
        }
    }

    if (minI == -1)
    {
        Info<< "Did not find face close to point " << nearPoint << endl;

        return -1;
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest face      : " << minI
        << " distance " << minDist
        << " to face centre " << mesh.faceCentres()[minI] << nl
        << "    next nearest face : " << almostMinI
        << " distance " << almostMinDist
        << " to face centre " << mesh.faceCentres()[almostMinI] << nl
        << endl;

    if (almostMinDist < 2*minDist)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}


// Find cell with cell centre close to given point.
label findCell(const primitiveMesh& mesh, const point& nearPoint)
{
    label celli = mesh.findCell(nearPoint);

    if (celli != -1)
    {
        scalar distToCcSqr = nearPoint.distSqr(mesh.cellCentres()[celli]);

        const labelList& cPoints = mesh.cellPoints()[celli];

        label minI = -1;
        scalar minDistSqr = GREAT;

        for (const label pointi : cPoints)
        {
            scalar distSqr = nearPoint.distSqr(mesh.points()[pointi]);

            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minI = pointi;
            }
        }

        // Decide if nearPoint unique enough.
        Info<< "Found to point " << nearPoint << nl
            << "    nearest cell       : " << celli
            << " distance " << Foam::sqrt(distToCcSqr)
            << " to cell centre " << mesh.cellCentres()[celli] << nl
            << "    nearest mesh point : " << minI
            << " distance " << Foam::sqrt(minDistSqr)
            << " to " << mesh.points()[minI] << nl
            << endl;

        if (minDistSqr < 4*distToCcSqr)
        {
            Info<< "Mesh point too close to nearest cell centre. Aborting"
                << endl;

            celli = -1;
        }
    }

    return celli;
}




int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Manipulate mesh elements.\n"
        "For example, moving points, splitting/collapsing edges etc."
    );
    #include "addOverwriteOption.H"
    argList::addOption("dict", "file", "Alternative modifyMeshDict");

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.found("overwrite");

    Info<< "Reading modifyMeshDict\n" << endl;

    // Read meshing dictionary
    const word dictName("modifyMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary dict(dictIO);

    // Read all from the dictionary.
    List<Pair<point>> pointsToMove(dict.lookup("pointsToMove"));
    List<Pair<point>> edgesToSplit(dict.lookup("edgesToSplit"));
    List<Pair<point>> facesToTriangulate
    (
        dict.lookup("facesToTriangulate")
    );
    // Optional
    List<Pair<point>> facesToSplit;
    dict.readIfPresent("facesToSplit", facesToSplit);

    bool cutBoundary =
    (
        pointsToMove.size()
     || edgesToSplit.size()
     || facesToTriangulate.size()
     || facesToSplit.size()
    );

    List<Pair<point>> edgesToCollapse(dict.lookup("edgesToCollapse"));

    bool collapseEdge = edgesToCollapse.size();

    List<Pair<point>> cellsToPyramidise(dict.lookup("cellsToSplit"));

    bool cellsToSplit = cellsToPyramidise.size();

    // List<Tuple2<pointField,point>>
    //     cellsToCreate(dict.lookup("cellsToCreate"));

    Info<< "Read from " << dict.name() << nl
        << "  Boundary cutting module:" << nl
        << "    points to move      :" << pointsToMove.size() << nl
        << "    edges to split      :" << edgesToSplit.size() << nl
        << "    faces to split      :" << facesToSplit.size() << nl
        << "    faces to triangulate:" << facesToTriangulate.size() << nl
        << "  Cell splitting module:" << nl
        << "    cells to split      :" << cellsToPyramidise.size() << nl
        << "  Edge collapsing module:" << nl
        << "    edges to collapse   :" << edgesToCollapse.size() << nl
        //<< "    cells to create     :" << cellsToCreate.size() << nl
        << endl;

    if
    (
        (cutBoundary && collapseEdge)
     || (cutBoundary && cellsToSplit)
     || (collapseEdge && cellsToSplit)
    )
    {
        FatalErrorInFunction
            << "Used more than one mesh modifying module "
            << "(boundary cutting, cell splitting, edge collapsing)" << nl
            << "Please do them in separate passes." << exit(FatalError);
    }



    // Get calculating engine for all of outside
    const SubList<face> outsideFaces
    (
        mesh.faces(),
        mesh.nBoundaryFaces(),
        mesh.nInternalFaces()
    );

    primitivePatch allBoundary(outsideFaces, mesh.points());


    // Look up mesh labels and convert to input for boundaryCutter.

    bool validInputs = true;


    Info<< nl << "Looking up points to move ..." << nl << endl;
    Map<point> pointToPos(pointsToMove.size());
    for (const Pair<point>& pts : pointsToMove)
    {
        const label pointi = findPoint(allBoundary, pts.first());

        if (pointi == -1 || !pointToPos.insert(pointi, pts.second()))
        {
            Info<< "Could not insert mesh point " << pointi
                << " for input point " << pts.first() << nl
                << "Perhaps the point is already marked for moving?" << endl;
            validInputs = false;
        }
    }


    Info<< nl << "Looking up edges to split ..." << nl << endl;
    Map<List<point>> edgeToCuts(edgesToSplit.size());
    for (const Pair<point>& pts : edgesToSplit)
    {
        label edgeI = findEdge(mesh, allBoundary, pts.first());

        if
        (
            edgeI == -1
        || !edgeToCuts.insert(edgeI, List<point>(1, pts.second()))
        )
        {
            Info<< "Could not insert mesh edge " << edgeI
                << " for input point " << pts.first() << nl
                << "Perhaps the edge is already marked for cutting?" << endl;

            validInputs = false;
        }
    }


    Info<< nl << "Looking up faces to triangulate ..." << nl << endl;
    Map<point> faceToDecompose(facesToTriangulate.size());
    for (const Pair<point>& pts : facesToTriangulate)
    {
        label facei = findFace(mesh, allBoundary, pts.first());

        if (facei == -1 || !faceToDecompose.insert(facei, pts.second()))
        {
            Info<< "Could not insert mesh face " << facei
                << " for input point " << pts.first() << nl
                << "Perhaps the face is already marked for splitting?" << endl;

            validInputs = false;
        }
    }


    Info<< nl << "Looking up faces to split ..." << nl << endl;
    Map<labelPair> faceToSplit(facesToSplit.size());
    for (const Pair<point>& pts : facesToSplit)
    {
        label facei = findFace(mesh, allBoundary, pts.first());
        if (facei == -1)
        {
            Info<< "Could not insert mesh face " << facei
                << " for input point " << pts.first() << nl
                << "Perhaps the face is already marked for splitting?" << endl;

            validInputs = false;
        }
        else
        {
            // Find nearest points on face
            const primitivePatch pp
            (
                SubList<face>(mesh.faces(), 1, facei),
                mesh.points()
            );

            const label p0 = findPoint(pp, pts.first());
            const label p1 = findPoint(pp, pts.second());

            const face& f = mesh.faces()[facei];

            if
            (
                p0 != -1
             && p1 != -1
             && faceToSplit.insert(facei, labelPair(f.find(p0), f.find(p1)))
            )
            {}
            else
            {
                Info<< "Could not insert mesh face " << facei
                    << " for input coordinates " << pts
                    << " with vertices " << p0 << " and " << p1 << nl
                    << "Perhaps the face is already marked for splitting?"
                    << endl;

            }
        }
    }


    Info<< nl << "Looking up cells to convert to pyramids around"
        << " cell centre ..." << nl << endl;
    Map<point> cellToPyrCentre(cellsToPyramidise.size());
    for (const Pair<point>& pts : cellsToPyramidise)
    {
        label celli = findCell(mesh, pts.first());

        if (celli == -1 || !cellToPyrCentre.insert(celli, pts.second()))
        {
            Info<< "Could not insert mesh cell " << celli
                << " for input point " << pts.first() << nl
                << "Perhaps the cell is already marked for splitting?" << endl;

            validInputs = false;
        }
    }


    Info<< nl << "Looking up edges to collapse ..." << nl << endl;
    Map<point> edgeToPos(edgesToCollapse.size());
    for (const Pair<point>& pts : edgesToCollapse)
    {
        label edgeI = findEdge(mesh, allBoundary, pts.first());

        if (edgeI == -1 || !edgeToPos.insert(edgeI, pts.second()))
        {
            Info<< "Could not insert mesh edge " << edgeI
                << " for input point " << pts.first() << nl
                << "Perhaps the edge is already marked for collaping?" << endl;

            validInputs = false;
        }
    }



    if (!validInputs)
    {
        Info<< nl << "There was a problem in one of the inputs in the"
            << " dictionary. Not modifying mesh." << endl;
    }
    else if (cellToPyrCentre.size())
    {
        Info<< nl << "All input cells located. Modifying mesh." << endl;

        // Mesh change engine
        cellSplitter cutter(mesh);

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Insert commands into meshMod
        cutter.setRefinement(cellToPyrCentre, meshMod);

        // Do changes
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        cutter.updateMesh(morphMap());

        if (!overwrite)
        {
            ++runTime;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Write resulting mesh
        Info<< "Writing modified mesh to time " << runTime.timeName() << endl;
        mesh.write();
        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }
    else if (edgeToPos.size())
    {
        Info<< nl << "All input edges located. Modifying mesh." << endl;

        // Mesh change engine
        edgeCollapser cutter(mesh);

        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        pointField newPoints(points);

        bitSet collapseEdge(mesh.nEdges());
        Map<point> collapsePointToLocation(mesh.nPoints());

        // Get new positions and construct collapse network
        forAllConstIters(edgeToPos, iter)
        {
            label edgeI = iter.key();
            const edge& e = edges[edgeI];

            collapseEdge.set(edgeI);
            collapsePointToLocation.set(e[1], points[e[0]]);

            newPoints[e[0]] = iter.val();
        }

        // Move master point to destination.
        mesh.movePoints(newPoints);

        List<pointEdgeCollapse> allPointInfo;
        const globalIndex globalPoints(mesh.nPoints());
        labelList pointPriority(mesh.nPoints(), Zero);

        cutter.consistentCollapse
        (
            globalPoints,
            pointPriority,
            collapsePointToLocation,
            collapseEdge,
            allPointInfo
        );

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Insert
        cutter.setRefinement(allPointInfo, meshMod);

        // Do changes
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Not implemented yet:
        //cutter.updateMesh(morphMap());


        if (!overwrite)
        {
            ++runTime;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Write resulting mesh
        Info<< "Writing modified mesh to time " << runTime.timeName() << endl;
        mesh.write();
        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }
    else
    {
        Info<< nl << "All input points located. Modifying mesh." << endl;

        // Mesh change engine
        boundaryCutter cutter(mesh);

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Insert commands into meshMod
        cutter.setRefinement
        (
            pointToPos,
            edgeToCuts,
            faceToSplit,        // Faces to split diagonally
            faceToDecompose,    // Faces to triangulate
            meshMod
        );

        // Do changes
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        cutter.updateMesh(morphMap());

        if (!overwrite)
        {
            ++runTime;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Write resulting mesh
        Info<< "Writing modified mesh to time " << runTime.timeName() << endl;
        mesh.write();
        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
