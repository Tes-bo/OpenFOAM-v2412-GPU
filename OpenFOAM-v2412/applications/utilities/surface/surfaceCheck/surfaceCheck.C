/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    surfaceCheck

Group
    grpSurfaceUtilities

Description
    Check geometric and topological quality of a surface.

Usage
    \b surfaceCheck [OPTION] surfaceFile

    Options:
      - \par -checkSelfIntersection
        Check for self-intersection.

      - \par -splitNonManifold
        Split surface along non-manifold edges.

      - \par -verbose
        Extra verbosity.

      - \par -blockMesh
        Write vertices/blocks for tight-fitting 1 cell blockMeshDict.

      - \par -outputThreshold \<num files\>
        Upper limit on the number of files written.
        Prevent surfaces with many disconnected parts from writing many files.
        The default is 10. A value of 0 suppresses file writing.

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "edgeHashes.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "SortableList.H"
#include "PatchTools.H"
#include "vtkSurfaceWriter.H"
#include "functionObject.H"
#include "DynamicField.H"
#include "edgeMesh.H"

using namespace Foam;

labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, Zero);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < SMALL)
        {
            index = 0;
        }
        else if (val >= max - SMALL)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((val - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                WarningInFunction
                    << "value " << val << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}



void writeZoning
(
    surfaceWriter& writer,
    const triSurface& surf,
    const labelList& faceZone,
    const word& fieldName,
    const fileName& surfFilePath,
    const word& surfFileStem
)
{
    // Transcribe faces
    faceList faces;
    surf.triFaceFaces(faces);

    writer.open
    (
        surf.points(),
        faces,
        (surfFilePath / surfFileStem),
        false  // serial - already merged
    );

    fileName outputName = writer.write(fieldName, labelField(faceZone));

    writer.clear();

    Info<< "Wrote zoning to " << outputName << nl << endl;
}


void writeParts
(
    const triSurface& surf,
    const label nFaceZones,
    const labelList& faceZone,
    const fileName& surfFilePath,
    const word& surfFileStem
)
{
    for (label zone = 0; zone < nFaceZones; zone++)
    {
        boolList includeMap(surf.size(), false);

        forAll(faceZone, facei)
        {
            if (faceZone[facei] == zone)
            {
                includeMap[facei] = true;
            }
        }

        triSurface subSurf(surf.subsetMesh(includeMap));

        fileName subName
        (
            surfFilePath
          / surfFileStem + "_" + name(zone) + ".obj"
        );

        Info<< "writing part " << zone << " size " << subSurf.size()
            << " to " << subName << endl;

        subSurf.write(subName);
    }
}


void syncEdges(const triSurface& p, labelHashSet& markedEdges)
{
    // See comment below about having duplicate edges

    const edgeList& edges = p.edges();
    edgeHashSet edgeSet(2*markedEdges.size());

    for (const label edgei : markedEdges)
    {
        edgeSet.insert(edges[edgei]);
    }

    forAll(edges, edgei)
    {
        if (edgeSet.found(edges[edgei]))
        {
            markedEdges.insert(edgei);
        }
    }
}


void syncEdges(const triSurface& p, boolList& isMarkedEdge)
{
    // See comment below about having duplicate edges

    const edgeList& edges = p.edges();

    label n = 0;
    forAll(isMarkedEdge, edgei)
    {
        if (isMarkedEdge[edgei])
        {
            n++;
        }
    }

    edgeHashSet edgeSet(2*n);

    forAll(isMarkedEdge, edgei)
    {
        if (isMarkedEdge[edgei])
        {
            edgeSet.insert(edges[edgei]);
        }
    }

    forAll(edges, edgei)
    {
        if (edgeSet.found(edges[edgei]))
        {
            isMarkedEdge[edgei] = true;
        }
    }
}


void writeEdgeSet
(
    const word& setName,
    const triSurface& surf,
    const labelUList& edgeSet
)
{
    // Get compact edge mesh
    labelList pointToCompact(surf.nPoints(), -1);
    DynamicField<point> compactPoints(edgeSet.size());
    DynamicList<edge> compactEdges(edgeSet.size());
    for (label edgei : edgeSet)
    {
        const edge& e = surf.edges()[edgei];
        edge compactEdge(-1, -1);
        forAll(e, ep)
        {
            label& compacti = pointToCompact[e[ep]];
            if (compacti == -1)
            {
                compacti = compactPoints.size();
                label pointi = surf.meshPoints()[e[ep]];
                compactPoints.append(surf.points()[pointi]);
            }
            compactEdge[ep] = compacti;
        }
        compactEdges.append(compactEdge);
    }

    edgeMesh eMesh(std::move(compactPoints), std::move(compactEdges));
    eMesh.write(setName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check geometric and topological quality of a surface"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");

    argList::addBoolOption
    (
        "checkSelfIntersection",
        "Also check for self-intersection"
    );
    argList::addBoolOption
    (
        "splitNonManifold",
        "Split surface along non-manifold edges "
        "(default split is fully disconnected)"
    );
    argList::addVerboseOption();
    argList::addBoolOption
    (
        "blockMesh",
        "Write vertices/blocks for blockMeshDict"
    );
    argList::addOption
    (
        "outputThreshold",
        "number",
        "Upper limit on the number of files written. "
        "Default is 10, using 0 suppresses file writing."
    );
    argList::addOption
    (
        "writeSets",
        "surfaceFormat",
        "Reconstruct and write problem triangles/edges in selected format"
    );


    argList args(argc, argv);

    const auto surfName = args.get<fileName>(1);
    const bool checkSelfIntersect = args.found("checkSelfIntersection");
    const bool splitNonManifold = args.found("splitNonManifold");
    const label outputThreshold =
        args.getOrDefault<label>("outputThreshold", 10);
    const word surfaceFormat = args.getOrDefault<word>("writeSets", "");
    const bool writeSets = !surfaceFormat.empty();

    autoPtr<surfaceWriter> surfWriter;
    word edgeFormat;
    if (writeSets)
    {
        surfWriter = surfaceWriter::New(surfaceFormat);

        // Option1: hard-coded format
        edgeFormat = "obj";
        //// Option2: same type as surface format. Problem is e.g. .obj format
        ////          is not a surfaceWriter format.
        //edgeFormat = surfaceFormat;
        //if (!edgeMesh::canWriteType(edgeFormat))
        //{
        //    edgeFormat = "obj";
        //}
    }


    Info<< "Reading surface from "
        << args.relativePath(surfName) << " ..." << nl << endl;

    // Read
    // ~~~~

    triSurface surf(surfName);


    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;


    // Split into path and stem (no extension)
    const fileName surfFilePath(surfName.path());
    word surfFileStem(surfName.stem());

    // If originally ".gz", need to strip extension again
    if (surfName.has_ext("gz"))
    {
        surfFileStem.remove_ext();
    }


    // write bounding box corners
    if (args.found("blockMesh"))
    {
        pointField cornerPts(boundBox(surf.points(), false).points());

        Info<< "// blockMeshDict info" << nl << nl;

        Info<< "vertices\n(" << nl;
        forAll(cornerPts, pti)
        {
            Info<< "    " << cornerPts[pti] << nl;
        }

        // number of divisions needs adjustment later
        Info<< ");\n" << nl
            << "blocks\n"
            << "(\n"
            << "    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)\n"
            << ");\n" << nl;

        Info<< "edges\n();" << nl
            << "patches\n();" << endl;

        Info<< nl << "// end blockMeshDict info" << nl << endl;
    }


    // Region sizes
    // ~~~~~~~~~~~~

    {
        labelList regionSize(surf.patches().size(), Foam::zero{});

        forAll(surf, facei)
        {
            label region = surf[facei].region();

            if (region < 0 || region >= regionSize.size())
            {
                WarningInFunction
                    << "Triangle " << facei << " vertices " << surf[facei]
                    << " has region " << region << " which is outside the range"
                    << " of regions 0.." << surf.patches().size()-1
                    << endl;
            }
            else
            {
                regionSize[region]++;
            }
        }

        Info<< "Index\tRegion\tSize" << nl
            << "-----\t------\t----" << nl;
        forAll(surf.patches(), patchi)
        {
            Info<< patchi << '\t'
                << surf.patches()[patchi].name() << '\t'
                << regionSize[patchi] << nl;
        }
        Info<< nl << endl;
    }


    // Check triangles
    // ~~~~~~~~~~~~~~~

    {
        DynamicList<label> illegalFaces(surf.size()/100 + 1);

        forAll(surf, facei)
        {
            // Check silently
            if (!triSurfaceTools::validTri(surf, facei, false))
            {
                illegalFaces.append(facei);
            }
        }

        if (illegalFaces.size())
        {
            Info<< "Surface has " << illegalFaces.size()
                << " illegal triangles." << endl;

            if (surfWriter)
            {
                boolList isIllegalFace(surf.size(), false);
                UIndirectList<bool>(isIllegalFace, illegalFaces) = true;

                triSurface subSurf(surf.subsetMesh(isIllegalFace));


                // Transcribe faces
                faceList faces;
                subSurf.triFaceFaces(faces);

                surfWriter->open
                (
                    subSurf.points(),
                    faces,
                    (surfFilePath / surfFileStem),
                    false // serial - already merged
                );

                fileName outputName = surfWriter->write
                (
                    "illegal",
                    scalarField(subSurf.size(), Zero)
                );

                surfWriter->clear();

                Info<< "Wrote illegal triangles to "
                    << outputName << nl << endl;
            }
            else if (outputThreshold > 0)
            {
                forAll(illegalFaces, i)
                {
                    // Force warning message
                    triSurfaceTools::validTri(surf, illegalFaces[i], true);
                    if (i >= outputThreshold)
                    {
                        Info<< "Suppressing further warning messages."
                            << " Use -outputThreshold to increase number"
                            << " of warnings." << endl;
                        break;
                    }
                }

                OFstream str("illegalFaces");
                Info<< "Dumping conflicting face labels to " << str.name()
                    << endl
                    << "Paste this into the input for surfaceSubset" << endl;
                str << illegalFaces;
            }
        }
        else
        {
            Info<< "Surface has no illegal triangles." << endl;
        }
        Info<< endl;
    }



    // Triangle quality
    // ~~~~~~~~~~~~~~~~

    {
        scalarField triQ(surf.size(), Zero);
        forAll(surf, facei)
        {
            const labelledTri& f = surf[facei];

            if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
            {
                //WarningInFunction
                //    << "Illegal triangle " << facei << " vertices " << f
                //    << " coords " << f.points(surf.points()) << endl;
            }
            else
            {
                triQ[facei] = triPointRef
                (
                    surf.points()[f[0]],
                    surf.points()[f[1]],
                    surf.points()[f[2]]
                ).quality();
            }
        }

        labelList binCount = countBins(0, 1, 20, triQ);

        Info<< "Triangle quality (equilateral=1, collapsed=0):"
            << endl;


        OSstream& os = Info;
        os.width(4);

        scalar dist = (1.0 - 0.0)/20.0;
        scalar min = 0;
        forAll(binCount, bini)
        {
            Info<< "    " << min << " .. " << min+dist << "  : "
                << 1.0/surf.size() * binCount[bini]
                << endl;
            min += dist;
        }
        Info<< endl;

        labelPair minMaxIds = findMinMax(triQ);

        const label minIndex = minMaxIds.first();
        const label maxIndex = minMaxIds.second();

        Info<< "    min " << triQ[minIndex] << " for triangle " << minIndex
            << nl
            << "    max " << triQ[maxIndex] << " for triangle " << maxIndex
            << nl
            << endl;


        if (triQ[minIndex] < SMALL)
        {
            WarningInFunction
                << "Minimum triangle quality is "
                << triQ[minIndex] << ". This might give problems in"
                << " self-intersection testing later on." << endl;
        }

        // Dump for subsetting
        if (surfWriter)
        {
            // Transcribe faces
            faceList faces(surf.size());
            surf.triFaceFaces(faces);

            surfWriter->open
            (
                surf.points(),
                faces,
                (surfFilePath / surfFileStem),
                false // serial - already merged
            );

            fileName outputName = surfWriter->write("quality", triQ);

            surfWriter->clear();

            Info<< "Wrote triangle-quality to "
                << outputName << nl << endl;
        }
        else if (outputThreshold > 0)
        {
            DynamicList<label> problemFaces(surf.size()/100+1);

            forAll(triQ, facei)
            {
                if (triQ[facei] < 1e-11)
                {
                    problemFaces.append(facei);
                }
            }

            if (!problemFaces.empty())
            {
                OFstream str("badFaces");

                Info<< "Dumping bad quality faces to " << str.name() << endl
                    << "Paste this into the input for surfaceSubset" << nl
                    << nl << endl;

                str << problemFaces;
            }
        }
    }



    // Edges
    // ~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        scalarField edgeMag(edges.size());

        forAll(edges, edgei)
        {
            edgeMag[edgei] = edges[edgei].mag(localPoints);
        }

        labelPair minMaxIds = findMinMax(edgeMag);

        const label minEdgei = minMaxIds.first();
        const label maxEdgei = minMaxIds.second();

        const edge& minE = edges[minEdgei];
        const edge& maxE = edges[maxEdgei];


        Info<< "Edges:" << nl
            << "    min " << edgeMag[minEdgei] << " for edge " << minEdgei
            << " points " << localPoints[minE[0]] << localPoints[minE[1]]
            << nl
            << "    max " << edgeMag[maxEdgei] << " for edge " << maxEdgei
            << " points " << localPoints[maxE[0]] << localPoints[maxE[1]]
            << nl
            << endl;
    }



    // Close points
    // ~~~~~~~~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        const boundBox bb(localPoints);
        scalar smallDim = 1e-6 * bb.mag();

        Info<< "Checking for points less than 1e-6 of bounding box ("
            << bb.span() << " metre) apart."
            << endl;

        // Sort points
        SortableList<scalar> sortedMag(mag(localPoints));

        label nClose = 0;

        for (label i = 1; i < sortedMag.size(); i++)
        {
            label pti = sortedMag.indices()[i];

            label prevPti = sortedMag.indices()[i-1];

            if (mag(localPoints[pti] - localPoints[prevPti]) < smallDim)
            {
                // Check if neighbours.
                const labelList& pEdges = surf.pointEdges()[pti];

                label edgei = -1;

                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if (e[0] == prevPti || e[1] == prevPti)
                    {
                        // point1 and point0 are connected through edge.
                        edgei = pEdges[i];

                        break;
                    }
                }

                if (nClose < outputThreshold)
                {
                    if (edgei == -1)
                    {
                        Info<< "    close unconnected points "
                            << pti << ' ' << localPoints[pti]
                            << " and " << prevPti << ' '
                            << localPoints[prevPti]
                            << " distance:"
                            << mag(localPoints[pti] - localPoints[prevPti])
                            << endl;
                    }
                    else
                    {
                        Info<< "    small edge between points "
                            << pti << ' ' << localPoints[pti]
                            << " and " << prevPti << ' '
                            << localPoints[prevPti]
                            << " distance:"
                            << mag(localPoints[pti] - localPoints[prevPti])
                            << endl;
                    }
                }
                else if (nClose == outputThreshold)
                {
                    Info<< "Suppressing further warning messages."
                        << " Use -outputThreshold to increase number"
                        << " of warnings." << endl;
                }

                nClose++;
            }
        }

        Info<< "Found " << nClose << " nearby points." << nl
            << endl;
    }



    // Check manifold
    // ~~~~~~~~~~~~~~

    DynamicList<label> problemFaces(surf.size()/100 + 1);
    DynamicList<label> openEdges(surf.nEdges()/100 + 1);
    DynamicList<label> multipleEdges(surf.nEdges()/100 + 1);

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgei)
    {
        const labelList& myFaces = edgeFaces[edgei];

        if (myFaces.size() == 1)
        {
            problemFaces.append(myFaces[0]);
            openEdges.append(edgei);
        }
    }

    forAll(edgeFaces, edgei)
    {
        const labelList& myFaces = edgeFaces[edgei];

        if (myFaces.size() > 2)
        {
            forAll(myFaces, myFacei)
            {
                problemFaces.append(myFaces[myFacei]);
            }
            multipleEdges.append(edgei);
        }
    }
    problemFaces.shrink();

    if (openEdges.size() || multipleEdges.size())
    {
        Info<< "Surface is not closed since not all edges ("
            << edgeFaces.size() << ") connected to "
            << "two faces:" << endl
            << "    connected to one face : " << openEdges.size() << endl
            << "    connected to >2 faces : " << multipleEdges.size() << endl;

        Info<< "Conflicting face labels:" << problemFaces.size() << endl;

        if (!edgeFormat.empty() && openEdges.size())
        {
            const fileName outputName
            (
                surfName.lessExt()
              + "_open."
              + edgeFormat
            );
            Info<< "Writing open edges to "
                << args.relativePath(outputName) << " ..." << endl;
            writeEdgeSet(outputName, surf, openEdges);
        }
        if (!edgeFormat.empty() && multipleEdges.size())
        {
            const fileName outputName
            (
                surfName.lessExt()
              + "_multiply."
              + edgeFormat
            );
            Info<< "Writing multiply connected edges to "
                << args.relativePath(outputName) << " ..." << endl;
            writeEdgeSet(outputName, surf, multipleEdges);
        }
    }
    else
    {
        Info<< "Surface is closed. All edges connected to two faces." << endl;
    }
    Info<< endl;



    // Check singly connected domain
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList borderEdge(surf.nEdges(), false);
        if (splitNonManifold)
        {
            forAll(edgeFaces, edgei)
            {
                if (edgeFaces[edgei].size() > 2)
                {
                    borderEdge[edgei] = true;
                }
            }
            syncEdges(surf, borderEdge);
        }

        labelList faceZone;
        label numZones = surf.markZones(borderEdge, faceZone);

        Info<< "Number of unconnected parts : " << numZones << endl;

        if (numZones > 1 && outputThreshold > 0)
        {
            Info<< "Splitting surface into parts ..." << endl << endl;

            if (!surfWriter)
            {
                surfWriter.reset(new surfaceWriters::vtkWriter());
            }

            writeZoning
            (
                *surfWriter,
                surf,
                faceZone,
                "zone",
                surfFilePath,
                surfFileStem
            );

            if (numZones > outputThreshold)
            {
                Info<< "Limiting number of files to " << outputThreshold
                    << endl;
            }
            writeParts
            (
                surf,
                min(outputThreshold, numZones),
                faceZone,
                surfFilePath,
                surfFileStem
            );
        }
    }



    // Check orientation
    // ~~~~~~~~~~~~~~~~~

    labelHashSet borderEdge(surf.size()/1000);
    PatchTools::checkOrientation(surf, false, &borderEdge);

    // Bit strange: if a triangle has two same vertices (illegal!) it will
    // still have three distinct edges (two of which have the same vertices).
    // In this case the faceEdges addressing is not symmetric, i.e. a
    // neighbouring, valid, triangle will have correct addressing so 3 distinct
    // edges so it will miss one of those two identical edges.
    // - we don't want to fix this in PrimitivePatch since it is too specific
    // - instead just make sure we mark all identical edges consistently
    //   when we use them for marking.

    syncEdges(surf, borderEdge);

    //
    // Colour all faces into zones using borderEdge
    //
    labelList normalZone;
    label numNormalZones = PatchTools::markZones(surf, borderEdge, normalZone);

    Info<< endl
        << "Number of zones (connected area with consistent normal) : "
        << numNormalZones << endl;

    if (numNormalZones > 1)
    {
        Info<< "More than one normal orientation." << endl;

        if (outputThreshold > 0)
        {
            if (!surfWriter)
            {
                surfWriter.reset(new surfaceWriters::vtkWriter());
            }

            writeZoning
            (
                *surfWriter,
                surf,
                normalZone,
                "normal",
                surfFilePath,
                surfFileStem
            );

            if (numNormalZones > outputThreshold)
            {
                Info<< "Limiting number of files to " << outputThreshold
                    << endl;
            }
            writeParts
            (
                surf,
                min(outputThreshold, numNormalZones),
                normalZone,
                surfFilePath,
                surfFileStem + "_normal"
            );
        }
    }
    Info<< endl;



    // Check self-intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (checkSelfIntersect)
    {
        Info<< "Checking self-intersection." << endl;

        triSurfaceSearch querySurf(surf);

        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

        autoPtr<OBJstream> intStreamPtr;
        if (outputThreshold > 0)
        {
            intStreamPtr.reset(new OBJstream("selfInterPoints.obj"));
        }

        label nInt = 0;

        forAll(surf.edges(), edgei)
        {
            const edge& e = surf.edges()[edgei];
            const point& start = surf.points()[surf.meshPoints()[e[0]]];
            const point& end = surf.points()[surf.meshPoints()[e[1]]];

            // Exclude hits of connected triangles
            treeDataTriSurface::findSelfIntersectOp exclOp(tree, edgei);

            pointIndexHit hitInfo(tree.findLineAny(start, end, exclOp));

            if (hitInfo.hit())
            {
                nInt++;

                if (intStreamPtr)
                {
                    intStreamPtr().write(hitInfo.point());
                }

                // Try and find from other side.
                pointIndexHit hitInfo2(tree.findLineAny(end, start, exclOp));

                if (hitInfo2.hit() && hitInfo.index() != hitInfo2.index())
                {
                    nInt++;

                    if (intStreamPtr)
                    {
                        intStreamPtr().write(hitInfo2.point());
                    }
                }
            }
        }

        //// Check very near triangles
        //{
        //    const pointField& localPoints = surf.localPoints();
        //
        //    const boundBox bb(localPoints);
        //    scalar smallDim = 1e-6 * bb.mag();
        //    scalar smallDimSqr = Foam::sqr(smallDim);
        //
        //    const pointField& faceCentres = surf.faceCentres();
        //    forAll(faceCentres, faceI)
        //    {
        //        const point& fc = faceCentres[faceI];
        //        pointIndexHit hitInfo
        //        (
        //            tree.findNearest
        //            (
        //                fc,
        //                smallDimSqr,
        //                findSelfNearOp(tree, faceI)
        //            )
        //        );
        //
        //        if (hitInfo.hit() && intStreamPtr)
        //        {
        //            intStreamPtr().write(hitInfo.point());
        //
        //            label nearFaceI = hitInfo.index();
        //            triPointRef nearTri(surf[nearFaceI].tri(surf.points()));
        //            triStreamPtr().write
        //            (
        //                surf[faceI].tri(surf.points()),
        //                false
        //            );
        //            triStreamPtr().write(nearTri, false);
        //            nInt++;
        //        }
        //    }
        //}


        if (nInt == 0)
        {
            Info<< "Surface is not self-intersecting" << endl;
        }
        else
        {
            Info<< "Surface is self-intersecting at " << nInt
                << " locations." << endl;

            if (intStreamPtr)
            {
                Info<< "Writing intersection points to "
                    << intStreamPtr().name() << endl;
            }
        }
        Info<< endl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
