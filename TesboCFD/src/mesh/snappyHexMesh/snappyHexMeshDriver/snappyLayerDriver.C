/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Description
    All to do with adding cell layers

\*----------------------------------------------------------------------------*/

#include "snappyLayerDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "meshRefinement.H"
#include "removePoints.H"
#include "pointFields.H"
#include "motionSmoother.H"
#include "unitConversion.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "addPatchCellLayer.H"
#include "mapDistributePolyMesh.H"
#include "OBJstream.H"
#include "layerParameters.H"
#include "combineFaces.H"
#include "IOmanip.H"
#include "globalIndex.H"
#include "DynamicField.H"
#include "PatchTools.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "zeroFixedValuePointPatchFields.H"
#include "calculatedPointPatchFields.H"
#include "cyclicSlipPointPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "localPointRegion.H"
#include "externalDisplacementMeshMover.H"
#include "scalarIOField.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappyLayerDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// For debugging: Dump displacement to .obj files
void Foam::snappyLayerDriver::dumpDisplacement
(
    const fileName& prefix,
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp,
    const List<extrudeMode>& extrudeStatus
)
{
    OBJstream dispStr(prefix + "_disp.obj");
    Info<< "Writing all displacements to " << dispStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        const point& pt = pp.localPoints()[patchPointi];
        dispStr.writeLine(pt, pt + patchDisp[patchPointi]);
    }


    OBJstream illStr(prefix + "_illegal.obj");
    Info<< "Writing invalid displacements to " << illStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        if (extrudeStatus[patchPointi] != EXTRUDE)
        {
            const point& pt = pp.localPoints()[patchPointi];
            illStr.writeLine(pt, pt + patchDisp[patchPointi]);
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::snappyLayerDriver::avgPointData
(
    const indirectPrimitivePatch& pp,
    const scalarField& pointFld
)
{
    auto tfaceFld = tmp<scalarField>::New(pp.size(), Zero);
    auto& faceFld = tfaceFld.ref();

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        if (f.size())
        {
            forAll(f, fp)
            {
                faceFld[facei] += pointFld[f[fp]];
            }
            faceFld[facei] /= f.size();
        }
    }
    return tfaceFld;
}


// Check that primitivePatch is not multiply connected. Collect non-manifold
// points in pointSet.
void Foam::snappyLayerDriver::checkManifold
(
    const indirectPrimitivePatch& fp,
    pointSet& nonManifoldPoints
)
{
    // Check for non-manifold points (surface pinched at point)
    fp.checkPointManifold(false, &nonManifoldPoints);

    // Check for edge-faces (surface pinched at edge)
    const labelListList& edgeFaces = fp.edgeFaces();

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            const edge& e = fp.edges()[edgei];

            nonManifoldPoints.insert(fp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(fp.meshPoints()[e[1]]);
        }
    }
}


void Foam::snappyLayerDriver::checkMeshManifold() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Checking mesh manifoldness ..." << endl;

    pointSet nonManifoldPoints
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints() / 100
    );

    // Build primitivePatch out of faces and check it for problems.
    checkManifold
    (
        indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                identity(mesh.boundaryMesh().range())  // All outside faces
            ),
            mesh.points()
        ),
        nonManifoldPoints
    );

    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nNonManif > 0)
    {
        Info<< "Outside of mesh is multiply connected across edges or"
            << " points." << nl
            << "This is not a fatal error but might cause some unexpected"
            << " behaviour." << nl
            //<< "Writing " << nNonManif
            //<< " points where this happens to pointSet "
            //<< nonManifoldPoints.name()
            << endl;

        //nonManifoldPoints.instance() = meshRefiner_.timeName();
        //nonManifoldPoints.write();
    }
    Info<< endl;
}



// Unset extrusion on point. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const label patchPointi,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointi] == EXTRUDE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }
    else if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }

    return false;
}


// Unset extrusion on face. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const face& localFace,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    bool unextruded = false;

    forAll(localFace, fp)
    {
        if
        (
            unmarkExtrusion
            (
                localFace[fp],
                patchDisp,
                patchNLayers,
                extrudeStatus
            )
        )
        {
            unextruded = true;
        }
    }
    return unextruded;
}


Foam::label Foam::snappyLayerDriver::constrainFp(const label sz, const label fp)
{
    if (fp >= sz)
    {
        return 0;
    }
    else if (fp < 0)
    {
        return sz-1;
    }
    else
    {
        return fp;
    }
}


void Foam::snappyLayerDriver::countCommonPoints
(
    const indirectPrimitivePatch& pp,
    const label facei,

    Map<label>& nCommonPoints
) const
{
    const faceList& localFaces = pp.localFaces();
    const labelListList& pointFaces = pp.pointFaces();

    const face& f = localFaces[facei];

    nCommonPoints.clear();

    forAll(f, fp)
    {
        label pointi = f[fp];
        const labelList& pFaces = pointFaces[pointi];

        forAll(pFaces, pFacei)
        {
            label nbFacei = pFaces[pFacei];

            if (facei < nbFacei)
            {
                // Only check once for each combination of two faces.
                ++(nCommonPoints(nbFacei, 0));
            }
        }
    }
}


bool Foam::snappyLayerDriver::checkCommonOrder
(
    const label nCommon,
    const face& curFace,
    const face& nbFace
) const
{
    forAll(curFace, fp)
    {
        // Get the index in the neighbouring face shared with curFace
        const label nb = nbFace.find(curFace[fp]);

        if (nb != -1)
        {

            // Check the whole face from nb onwards for shared vertices
            // with neighbouring face. Rule is that any shared vertices
            // should be consecutive on both faces i.e. if they are
            // vertices fp,fp+1,fp+2 on one face they should be
            // vertices nb, nb+1, nb+2 (or nb+2, nb+1, nb) on the
            // other face.


            // Vertices before and after on curFace
            label fpPlus1 = curFace.fcIndex(fp);
            label fpMin1  = curFace.rcIndex(fp);

            // Vertices before and after on nbFace
            label nbPlus1 = nbFace.fcIndex(nb);
            label nbMin1  = nbFace.rcIndex(nb);

            // Find order of walking by comparing next points on both
            // faces.
            label curInc = labelMax;
            label nbInc = labelMax;

            if (nbFace[nbPlus1] == curFace[fpPlus1])
            {
                curInc = 1;
                nbInc = 1;
            }
            else if (nbFace[nbPlus1] == curFace[fpMin1])
            {
                curInc = -1;
                nbInc = 1;
            }
            else if (nbFace[nbMin1] == curFace[fpMin1])
            {
                curInc = -1;
                nbInc = -1;
            }
            else
            {
                curInc = 1;
                nbInc = -1;
            }


            // Pass1: loop until start of common vertices found.
            label curNb = nb;
            label curFp = fp;

            do
            {
                curFp = constrainFp(curFace.size(), curFp+curInc);
                curNb = constrainFp(nbFace.size(), curNb+nbInc);
            } while (curFace[curFp] == nbFace[curNb]);

            // Pass2: check equality walking from curFp, curNb
            // in opposite order.

            curInc = -curInc;
            nbInc = -nbInc;

            for (label commonI = 0; commonI < nCommon; commonI++)
            {
                curFp = constrainFp(curFace.size(), curFp+curInc);
                curNb = constrainFp(nbFace.size(), curNb+nbInc);

                if (curFace[curFp] != nbFace[curNb])
                {
                    // Error: gap in string of connected vertices
                    return false;
                }
            }

            // Done the curFace - nbFace combination.
            break;
        }
    }

    return true;
}


void Foam::snappyLayerDriver::checkCommonOrder
(
    const indirectPrimitivePatch& pp,
    const label facei,
    const Map<label>& nCommonPoints,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    forAllConstIters(nCommonPoints, iter)
    {
        const label nbFacei = iter.key();
        const label nCommon = iter.val();

        const face& curFace = pp[facei];
        const face& nbFace = pp[nbFacei];

        if
        (
            nCommon >= 2
         && nCommon != nbFace.size()
         && nCommon != curFace.size()
        )
        {
            bool stringOk = checkCommonOrder(nCommon, curFace, nbFace);

            if (!stringOk)
            {
                // Note: unmark whole face or just the common points?
                // For now unmark the whole face
                unmarkExtrusion
                (
                    pp.localFaces()[facei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
                unmarkExtrusion
                (
                    pp.localFaces()[nbFacei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }
}


void Foam::snappyLayerDriver::handleNonStringConnected
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    // Detect faces which are connected on non-consecutive vertices.
    // This is the "<<Number of faces with non-consecutive shared points"
    // warning from checkMesh. These faces cannot be extruded so
    // there is no need to even attempt it.

    List<extrudeMode> oldExtrudeStatus;
    autoPtr<OBJstream> str;
    if (debug&meshRefinement::LAYERINFO)
    {
        oldExtrudeStatus = extrudeStatus;
        str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
               /"nonStringConnected.obj"
            )
        );
        Pout<< "Dumping string edges to " << str().name();
    }


    // 1) Local
    Map<label> nCommonPoints(128);

    forAll(pp, facei)
    {
        countCommonPoints(pp, facei, nCommonPoints);

        // Faces share pointi. Find any more shared points
        // and if not in single string unmark all. See
        // primitiveMesh::checkCommonOrder
        checkCommonOrder
        (
            pp,
            facei,
            nCommonPoints,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );
    }

    // 2) TDB. Other face remote



    if (debug&meshRefinement::LAYERINFO)
    {
        forAll(extrudeStatus, pointi)
        {
            if (extrudeStatus[pointi] != oldExtrudeStatus[pointi])
            {
                str().write
                (
                    meshRefiner_.mesh().points()[pp.meshPoints()[pointi]]
                );
            }
        }
    }
}


// No extrusion at non-manifold points.
void Foam::snappyLayerDriver::handleNonManifolds
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const labelListList& edgeGlobalFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling non-manifold points ..." << endl;

    // Detect non-manifold points
    Info<< nl << "Checking patch manifoldness ..." << endl;

    pointSet nonManifoldPoints(mesh, "nonManifoldPoints", pp.nPoints());

    // 1. Local check. Note that we do not check for e.g. two patch faces
    // being connected via a point since their connection might be ok
    // through a coupled patch. The ultimate is to do a proper point-face
    // walk which is done when actually duplicating the points. Here we just
    // do the obvious problems.
    {
        // Check for edge-faces (surface pinched at edge)
        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = edgeFaces[edgei];
            if (eFaces.size() > 2)
            {
                const edge& e = pp.edges()[edgei];
                nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
                nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
            }
        }
    }

    // 2. Remote check for boundary edges on coupled boundaries
    forAll(edgeGlobalFaces, edgei)
    {
        if (edgeGlobalFaces[edgei].size() > 2)
        {
            // So boundary edges that are connected to more than 2 processors
            // i.e. a non-manifold edge which is exactly on a processor
            // boundary.
            const edge& e = pp.edges()[edgei];
            nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
        }
    }


    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    Info<< "Outside of local patch is multiply connected across edges or"
        << " points at " << nNonManif << " points." << endl;

    if (nNonManif > 0)
    {
        // Make sure all processors use the same information. The edge might
        // not exist locally but remotely there might be a problem with this
        // edge.
        nonManifoldPoints.sync(mesh);

        const labelList& meshPoints = pp.meshPoints();

        forAll(meshPoints, patchPointi)
        {
            if (nonManifoldPoints.found(meshPoints[patchPointi]))
            {
                unmarkExtrusion
                (
                    patchPointi,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }

    Info<< "Set displacement to zero for all " << nNonManif
        << " non-manifold points" << endl;



    // 4. Check for extrusion of baffles i.e. all edges of a face having the
    //    same two neighbouring faces (one of which is the current face).
    //    Note: this is detected locally already before - this test is for the
    //          extremely rare occurrence where the baffle faces are on
    //          different processors.
    {
        label nBaffleFaces = 0;

        const labelListList& faceEdges = pp.faceEdges();
        forAll(pp, facei)
        {
            const labelList& fEdges = faceEdges[facei];

            const labelList& globFaces0 = edgeGlobalFaces[fEdges[0]];
            if (globFaces0.size() == 2)
            {
                const edge e0(globFaces0[0], globFaces0[1]);
                bool isBaffle = true;
                for (label fp = 1; fp < fEdges.size(); fp++)
                {
                    const labelList& globFaces = edgeGlobalFaces[fEdges[fp]];
                    if
                    (
                        (globFaces.size() != 2)
                     || (edge(globFaces[0], globFaces[1]) != e0)
                    )
                    {
                        isBaffle = false;
                        break;
                    }
                }

                if (isBaffle)
                {
                    bool unextrude = unmarkExtrusion
                    (
                        pp.localFaces()[facei],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    if (unextrude)
                    {
                        //Pout<< "Detected extrusion of baffle face "
                        //    << pp.faceCentres()[facei]
                        //    << " since all edges have the same neighbours "
                        //    << e0 << endl;

                        nBaffleFaces++;
                    }
                }
            }
        }

        reduce(nBaffleFaces, sumOp<label>());

        if (nBaffleFaces)
        {
            Info<< "Set displacement to zero for all points on " << nBaffleFaces
                << " baffle faces" << endl;
        }
    }
}


// Parallel feature edge detection. Assumes non-manifold edges already handled.
void Foam::snappyLayerDriver::handleFeatureAngle
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const scalar minAngle,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const scalar minCos = Foam::cos(degToRad(minAngle));

    Info<< nl << "Handling feature edges (angle < " << minAngle
        << ") ..." << endl;

    if (minCos < 1-SMALL && minCos > -1+SMALL)
    {
        // Normal component of normals of connected faces.
        vectorField edgeNormal(mesh.nEdges(), point::max);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = pp.edgeFaces()[edgei];

            label meshEdgei = meshEdges[edgei];

            forAll(eFaces, i)
            {
                nomalsCombine()
                (
                    edgeNormal[meshEdgei],
                    pp.faceNormals()[eFaces[i]]
                );
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeNormal,
            nomalsCombine(),
            point::max          // null value
        );

        autoPtr<OBJstream> str;
        if (debug&meshRefinement::MESH)
        {
            str.reset
            (
                new OBJstream
                (
                    mesh.time().path()
                  / "featureEdges_"
                  + meshRefiner_.timeName()
                  + ".obj"
                )
            );
            Info<< "Writing feature edges to " << str().name() << endl;
        }

        label nFeats = 0;

        // Now on coupled edges the edgeNormal will have been truncated and
        // only be still be the old value where two faces have the same normal
        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = pp.edgeFaces()[edgei];

            label meshEdgei = meshEdges[edgei];

            const vector& n = edgeNormal[meshEdgei];

            if (n != point::max)
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];

                if (cos < minCos)
                {
                    const edge& e = pp.edges()[edgei];

                    unmarkExtrusion
                    (
                        e[0],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    unmarkExtrusion
                    (
                        e[1],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );

                    nFeats++;

                    if (str)
                    {
                        str().write(e, pp.localPoints());
                    }
                }
            }
        }

        Info<< "Set displacement to zero for points on "
            << returnReduce(nFeats, sumOp<label>())
            << " feature edges" << endl;
    }
}


// No extrusion on cells with warped faces. Calculates the thickness of the
// layer and compares it to the space the warped face takes up. Disables
// extrusion if layer thickness is more than faceRatio of the thickness of
// the face.
void Foam::snappyLayerDriver::handleWarpedFaces
(
    const indirectPrimitivePatch& pp,
    const scalar faceRatio,
    const boolList& relativeSizes,
    const scalar edge0Len,
    const labelList& cellLevel,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    Info<< nl << "Handling cells with warped patch faces ..." << nl;

    const pointField& points = mesh.points();

    // Local reference to face centres also used to trigger consistent
    // [re-]building of demand-driven face centres and areas
    const vectorField& faceCentres = mesh.faceCentres();

    label nWarpedFaces = 0;

    forAll(pp, i)
    {
        const face& f = pp[i];
        label faceI = pp.addressing()[i];
        label patchI = patches.patchID(faceI);

        // It is hard to calculate some length scale if not in relative
        // mode so disable this check.

        if (relativeSizes[patchI] && f.size() > 3)
        {
            label ownLevel = cellLevel[mesh.faceOwner()[faceI]];
            scalar edgeLen = edge0Len/(1<<ownLevel);

            // Normal distance to face centre plane
            const point& fc = faceCentres[faceI];
            const vector& fn = pp.faceNormals()[i];

            scalarField vProj(f.size());

            forAll(f, fp)
            {
                vector n = points[f[fp]] - fc;
                vProj[fp] = (n & fn);
            }

            // Get normal 'span' of face
            scalar minVal = min(vProj);
            scalar maxVal = max(vProj);

            if ((maxVal - minVal) > faceRatio * edgeLen)
            {
                if
                (
                    unmarkExtrusion
                    (
                        pp.localFaces()[i],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nWarpedFaces++;
                }
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nWarpedFaces, sumOp<label>())
        << " warped faces since layer would be > " << faceRatio
        << " of the size of the bounding box." << endl;
}


//// No extrusion on cells with multiple patch faces. There usually is a reason
//// why combinePatchFaces hasn't succeeded.
//void Foam::snappyLayerDriver::handleMultiplePatchFaces
//(
//    const indirectPrimitivePatch& pp,
//    pointField& patchDisp,
//    labelList& patchNLayers,
//    List<extrudeMode>& extrudeStatus
//) const
//{
//    const fvMesh& mesh = meshRefiner_.mesh();
//
//    Info<< nl << "Handling cells with multiple patch faces ..." << nl;
//
//    const labelListList& pointFaces = pp.pointFaces();
//
//    // Cells that should not get an extrusion layer
//    cellSet multiPatchCells(mesh, "multiPatchCells", pp.size());
//
//    // Detect points that use multiple faces on same cell.
//    forAll(pointFaces, patchPointi)
//    {
//        const labelList& pFaces = pointFaces[patchPointi];
//
//        labelHashSet pointCells(pFaces.size());
//
//        forAll(pFaces, i)
//        {
//            label celli = mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//            if (!pointCells.insert(celli))
//            {
//                // Second or more occurrence of cell so cell has two or more
//                // pp faces connected to this point.
//                multiPatchCells.insert(celli);
//            }
//        }
//    }
//
//    label nMultiPatchCells = returnReduce
//    (
//        multiPatchCells.size(),
//        sumOp<label>()
//    );
//
//    Info<< "Detected " << nMultiPatchCells
//        << " cells with multiple (connected) patch faces." << endl;
//
//    label nChanged = 0;
//
//    if (nMultiPatchCells > 0)
//    {
//        multiPatchCells.instance() = meshRefiner_.timeName();
//        Info<< "Writing " << nMultiPatchCells
//            << " cells with multiple (connected) patch faces to cellSet "
//            << multiPatchCells.objectPath() << endl;
//        multiPatchCells.write();
//
//
//        // Go through all points and remove extrusion on any cell in
//        // multiPatchCells
//        // (has to be done in separate loop since having one point on
//        // multipatches has to reset extrusion on all points of cell)
//
//        forAll(pointFaces, patchPointi)
//        {
//            if (extrudeStatus[patchPointi] != NOEXTRUDE)
//            {
//                const labelList& pFaces = pointFaces[patchPointi];
//
//                forAll(pFaces, i)
//                {
//                    label celli =
//                        mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//                    if (multiPatchCells.found(celli))
//                    {
//                        if
//                        (
//                            unmarkExtrusion
//                            (
//                                patchPointi,
//                                patchDisp,
//                                patchNLayers,
//                                extrudeStatus
//                            )
//                        )
//                        {
//                            nChanged++;
//                        }
//                    }
//                }
//            }
//        }
//
//        reduce(nChanged, sumOp<label>());
//    }
//
//    Info<< "Prevented extrusion on " << nChanged
//        << " points due to multiple patch faces." << nl << endl;
//}


void Foam::snappyLayerDriver::setNumLayers
(
    meshRefinement& meshRefiner,
    const labelList& patchToNLayers,
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    label& nAddedCells
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    Info<< nl << "Handling points with inconsistent layer specification ..."
        << endl;

    // Get for every point (really only necessary on patch external points)
    // the max and min of any patch faces using it.
    labelList maxLayers(patchNLayers.size(), labelMin);
    labelList minLayers(patchNLayers.size(), labelMax);

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        const labelList& meshPoints = mesh.boundaryMesh()[patchi].meshPoints();

        label wantedLayers = patchToNLayers[patchi];

        forAll(meshPoints, patchPointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[patchPointi]];

            maxLayers[ppPointi] = max(wantedLayers, maxLayers[ppPointi]);
            minLayers[ppPointi] = min(wantedLayers, minLayers[ppPointi]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxLayers,
        maxEqOp<label>(),
        labelMin            // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minLayers,
        minEqOp<label>(),
        labelMax            // null value
    );

    // Unmark any point with different min and max
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //label nConflicts = 0;

    forAll(maxLayers, i)
    {
        if (maxLayers[i] == labelMin || minLayers[i] == labelMax)
        {
            FatalErrorInFunction
                << "Patchpoint:" << i << " coord:" << pp.localPoints()[i]
                << " maxLayers:" << maxLayers
                << " minLayers:" << minLayers
                << abort(FatalError);
        }
        else if (maxLayers[i] == minLayers[i])
        {
            // Ok setting.
            patchNLayers[i] = maxLayers[i];
        }
        else
        {
            // Inconsistent num layers between patch faces using point
            //if
            //(
            //    unmarkExtrusion
            //    (
            //        i,
            //        patchDisp,
            //        patchNLayers,
            //        extrudeStatus
            //    )
            //)
            //{
            //    nConflicts++;
            //}
            patchNLayers[i] = maxLayers[i];
        }
    }


    // Calculate number of cells to create
    nAddedCells = 0;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        // Get max of extrusion per point
        label nCells = 0;
        forAll(f, fp)
        {
            nCells = max(nCells, patchNLayers[f[fp]]);
        }

        nAddedCells += nCells;
    }
    reduce(nAddedCells, sumOp<label>());

    //reduce(nConflicts, sumOp<label>());
    //
    //Info<< "Set displacement to zero for " << nConflicts
    //    << " points due to points being on multiple regions"
    //    << " with inconsistent nLayers specification." << endl;
}


// Construct pointVectorField with correct boundary conditions for adding
// layers
Foam::tmp<Foam::pointVectorField>
Foam::snappyLayerDriver::makeLayerDisplacementField
(
    const pointMesh& pMesh,
    const labelList& numLayers
)
{
    // Construct displacement field.
    const pointBoundaryMesh& pointPatches = pMesh.boundary();

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchVectorField::typeName
    );
    wordList actualPatchTypes(patchFieldTypes.size());
    forAll(pointPatches, patchi)
    {
        actualPatchTypes[patchi] = pointPatches[patchi].type();
    }

    forAll(numLayers, patchi)
    {
        //  0 layers: do not allow slip so fixedValue 0
        // >0 layers: fixedValue which gets adapted
        if (numLayers[patchi] == 0)
        {
            patchFieldTypes[patchi] =
                zeroFixedValuePointPatchVectorField::typeName;
        }
        else if (numLayers[patchi] > 0)
        {
            patchFieldTypes[patchi] = fixedValuePointPatchVectorField::typeName;
        }
    }

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
        }
    }


    const polyMesh& mesh = pMesh();

    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.

    return tmp<pointVectorField>::New
    (
        IOobject
        (
            "pointDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector(dimLength, Zero),
        patchFieldTypes,
        actualPatchTypes
    );
}


void Foam::snappyLayerDriver::growNoExtrusion
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Growing non-extrusion points by one layer ..." << endl;

    List<extrudeMode> grownExtrudeStatus(extrudeStatus);

    const faceList& localFaces = pp.localFaces();

    label nGrown = 0;

    forAll(localFaces, facei)
    {
        const face& f = localFaces[facei];

        bool hasSqueeze = false;
        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                hasSqueeze = true;
                break;
            }
        }

        if (hasSqueeze)
        {
            // Squeeze all points of face
            forAll(f, fp)
            {
                if
                (
                    extrudeStatus[f[fp]] == EXTRUDE
                 && grownExtrudeStatus[f[fp]] != NOEXTRUDE
                )
                {
                    grownExtrudeStatus[f[fp]] = NOEXTRUDE;
                    nGrown++;
                }
            }
        }
    }

    extrudeStatus.transfer(grownExtrudeStatus);


    // Synchronise since might get called multiple times.
    // Use the fact that NOEXTRUDE is the minimum value.
    {
        labelList status(extrudeStatus.size());
        forAll(status, i)
        {
            status[i] = extrudeStatus[i];
        }
        syncTools::syncPointList
        (
            meshRefiner_.mesh(),
            pp.meshPoints(),
            status,
            minEqOp<label>(),
            labelMax            // null value
        );
        forAll(status, i)
        {
            extrudeStatus[i] = extrudeMode(status[i]);
        }
    }


    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }

    reduce(nGrown, sumOp<label>());

    Info<< "Set displacement to zero for an additional " << nGrown
        << " points." << endl;
}


void Foam::snappyLayerDriver::determineSidePatches
(
    meshRefinement& meshRefiner,
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,

    labelList& edgePatchID,
    labelList& edgeZoneID,
    boolList& edgeFlip,
    labelList& inflateFaceID
)
{
    // Sometimes edges-to-be-extruded are on more than 2 processors.
    // Work out which 2 hold the faces to be extruded and thus which procpatch
    // the edge-face should be in. As an additional complication this might
    // mean that 2 procesors that were only edge-connected now suddenly need
    // to become face-connected i.e. have a processor patch between them.

    fvMesh& mesh = meshRefiner.mesh();

    // Determine edgePatchID. Any additional processor boundary gets added to
    // patchToNbrProc,nbrProcToPatch and nPatches gets set to the new number
    // of patches.
    // Note: faceZones are at this point split into baffles so any zone
    //       information might also come from boundary faces (hence
    //       zoneFromAnyFace set in call to calcExtrudeInfo)
    label nPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;
    addPatchCellLayer::calcExtrudeInfo
    (
        true,           // zoneFromAnyFace

        mesh,
        globalFaces,
        edgeGlobalFaces,
        pp,

        edgePatchID,
        nPatches,
        nbrProcToPatch,
        patchToNbrProc,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );

    label nOldPatches = mesh.boundaryMesh().size();
    label nAdded = returnReduce(nPatches-nOldPatches, sumOp<label>());
    Info<< nl << "Adding in total " << nAdded/2 << " inter-processor patches to"
        << " handle extrusion of non-manifold processor boundaries."
        << endl;

    if (nAdded > 0)
    {
        // We might not add patches in same order as in patchToNbrProc
        // so prepare to renumber edgePatchID
        Map<label> wantedToAddedPatch;

        for (label patchi = nOldPatches; patchi < nPatches; patchi++)
        {
            label nbrProci = patchToNbrProc[patchi];
            word name
            (
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProci)
            );

            dictionary patchDict;
            patchDict.add("type", processorPolyPatch::typeName);
            patchDict.add("myProcNo", Pstream::myProcNo());
            patchDict.add("neighbProcNo", nbrProci);
            patchDict.add("nFaces", 0);
            patchDict.add("startFace", mesh.nFaces());

            //Pout<< "Adding patch " << patchi
            //    << " name:" << name
            //    << " between " << Pstream::myProcNo()
            //    << " and " << nbrProci << endl;

            label procPatchi = meshRefiner.appendPatch
            (
                mesh,
                mesh.boundaryMesh().size(), // new patch index
                name,
                patchDict
            );
            wantedToAddedPatch.insert(patchi, procPatchi);
        }

        // Renumber edgePatchID
        forAll(edgePatchID, i)
        {
            label patchi = edgePatchID[i];
            const auto fnd = wantedToAddedPatch.cfind(patchi);
            if (fnd.good())
            {
                edgePatchID[i] = fnd.val();
            }
        }

        mesh.clearOut();
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
    }
}


void Foam::snappyLayerDriver::calculateLayerThickness
(
    const indirectPrimitivePatch& pp,
    const labelList& patchIDs,
    const layerParameters& layerParams,
    const labelList& cellLevel,
    const labelList& patchNLayers,
    const scalar edge0Len,

    scalarField& thickness,
    scalarField& minThickness,
    scalarField& expansionRatio
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Rework patch-wise layer parameters into minimum per point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: only layer parameters consistent with layer specification
    // method (see layerParameters) will be correct.
    scalarField firstLayerThickness(pp.nPoints(), GREAT);
    scalarField finalLayerThickness(pp.nPoints(), GREAT);
    scalarField totalThickness(pp.nPoints(), GREAT);
    scalarField expRatio(pp.nPoints(), GREAT);

    minThickness.setSize(pp.nPoints());
    minThickness = GREAT;

    thickness.setSize(pp.nPoints());
    thickness = GREAT;

    expansionRatio.setSize(pp.nPoints());
    expansionRatio = GREAT;

    for (const label patchi : patchIDs)
    {
        const labelList& meshPoints = patches[patchi].meshPoints();

        forAll(meshPoints, patchPointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[patchPointi]];

            firstLayerThickness[ppPointi] = min
            (
                firstLayerThickness[ppPointi],
                layerParams.firstLayerThickness()[patchi]
            );
            finalLayerThickness[ppPointi] = min
            (
                finalLayerThickness[ppPointi],
                layerParams.finalLayerThickness()[patchi]
            );
            totalThickness[ppPointi] = min
            (
                totalThickness[ppPointi],
                layerParams.thickness()[patchi]
            );
            expRatio[ppPointi] = min
            (
                expRatio[ppPointi],
                layerParams.expansionRatio()[patchi]
            );
            minThickness[ppPointi] = min
            (
                minThickness[ppPointi],
                layerParams.minThickness()[patchi]
            );
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        firstLayerThickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        finalLayerThickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        totalThickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        expRatio,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minThickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );


    // Now the thicknesses are set according to the minimum of connected
    // patches.

    // Determine per point the max cell level of connected cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList maxPointLevel(pp.nPoints(), labelMin);
    {
        forAll(pp, i)
        {
            label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];

            const face& f = pp.localFaces()[i];

            forAll(f, fp)
            {
                maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );
    }


    // Rework relative thickness into absolute
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // by multiplying with the internal cell size.
    // Note that we cannot loop over the patches since then points on
    // multiple patches would get multiplied with edgeLen twice ..
    {
        // Multiplication factor for relative sizes
        scalarField edgeLen(pp.nPoints(), GREAT);

        labelList spec(pp.nPoints(), layerParameters::FIRST_AND_TOTAL);

        bitSet isRelativePoint(mesh.nPoints());

        for (const label patchi : patchIDs)
        {
            const labelList& meshPoints = patches[patchi].meshPoints();
            const layerParameters::thicknessModelType patchSpec =
                layerParams.layerModels()[patchi];
            const bool relSize = layerParams.relativeSizes()[patchi];

            for (const label meshPointi : meshPoints)
            {
                const label ppPointi = pp.meshPointMap()[meshPointi];

                // Note: who wins if different specs?

                // Calculate undistorted edge size for this level.
                edgeLen[ppPointi] = min
                (
                    edgeLen[ppPointi],
                    edge0Len/(1<<maxPointLevel[ppPointi])
                );
                spec[ppPointi] = max(spec[ppPointi], patchSpec);
                isRelativePoint[meshPointi] =
                    isRelativePoint[meshPointi]
                 || relSize;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            edgeLen,
            minEqOp<scalar>(),
            GREAT                               // null value
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            spec,
            maxEqOp<label>(),
            label(layerParameters::FIRST_AND_TOTAL) // null value
        );
        syncTools::syncPointList
        (
            mesh,
            isRelativePoint,
            orEqOp<unsigned int>(),
            0
        );




        forAll(pp.meshPoints(), pointi)
        {
            const label meshPointi = pp.meshPoints()[pointi];
            const layerParameters::thicknessModelType pointSpec =
                static_cast<layerParameters::thicknessModelType>(spec[pointi]);

            if (pointSpec == layerParameters::FIRST_AND_RELATIVE_FINAL)
            {
                // This overrules the relative sizes flag for
                // first (always absolute) and final (always relative)
                finalLayerThickness[pointi] *= edgeLen[pointi];
                if (isRelativePoint[meshPointi])
                {
                    totalThickness[pointi] *= edgeLen[pointi];
                    minThickness[pointi] *= edgeLen[pointi];
                }
            }
            else if (isRelativePoint[meshPointi])
            {
                firstLayerThickness[pointi] *= edgeLen[pointi];
                finalLayerThickness[pointi] *= edgeLen[pointi];
                totalThickness[pointi] *= edgeLen[pointi];
                minThickness[pointi] *= edgeLen[pointi];
            }

            thickness[pointi] = min
            (
                thickness[pointi],
                layerParameters::layerThickness
                (
                    pointSpec,
                    patchNLayers[pointi],
                    firstLayerThickness[pointi],
                    finalLayerThickness[pointi],
                    totalThickness[pointi],
                    expRatio[pointi]
                )
            );
            expansionRatio[pointi] = min
            (
                expansionRatio[pointi],
                layerParameters::layerExpansionRatio
                (
                    pointSpec,
                    patchNLayers[pointi],
                    firstLayerThickness[pointi],
                    finalLayerThickness[pointi],
                    totalThickness[pointi],
                    expRatio[pointi]
                )
            );
        }
    }

    // Synchronise the determined thicknes. Note that this should not be
    // necessary since the inputs to the calls to layerThickness,
    // layerExpansionRatio above are already parallel consistent

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        thickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        expansionRatio,
        minEqOp<scalar>(),
        GREAT               // null value
    );

    //Info<< "calculateLayerThickness : min:" << gMin(thickness)
    //    << " max:" << gMax(thickness) << endl;

    // Print a bit
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        const int oldPrecision = Info.stream().precision();

        // Find maximum length of a patch name, for a nicer output
        label maxPatchNameLen = 0;
        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];
            word patchName = patches[patchi].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }

        Info<< nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
            << setw(0) << " faces    layers avg thickness[m]" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << " "
            << setw(0) << "                 near-wall overall" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
            << setw(0) << " -----    ------ --------- -------" << endl;


        const bitSet isMasterPoint(syncTools::getMasterPoints(mesh));

        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];

            const labelList& meshPoints = patches[patchi].meshPoints();
            const layerParameters::thicknessModelType spec =
                layerParams.layerModels()[patchi];

            scalar sumThickness = 0;
            scalar sumNearWallThickness = 0;
            label nMasterPoints = 0;

            forAll(meshPoints, patchPointi)
            {
                label meshPointi = meshPoints[patchPointi];
                if (isMasterPoint[meshPointi])
                {
                    label ppPointi = pp.meshPointMap()[meshPointi];

                    sumThickness += thickness[ppPointi];
                    sumNearWallThickness += layerParams.firstLayerThickness
                    (
                        spec,
                        patchNLayers[ppPointi],
                        firstLayerThickness[ppPointi],
                        finalLayerThickness[ppPointi],
                        thickness[ppPointi],
                        expansionRatio[ppPointi]
                    );
                    nMasterPoints++;
                }
            }

            label totNPoints = returnReduce(nMasterPoints, sumOp<label>());

            // For empty patches, totNPoints is 0.
            scalar avgThickness = 0;
            scalar avgNearWallThickness = 0;

            if (totNPoints > 0)
            {
                avgThickness =
                    returnReduce(sumThickness, sumOp<scalar>())
                  / totNPoints;
                avgNearWallThickness =
                    returnReduce(sumNearWallThickness, sumOp<scalar>())
                  / totNPoints;
            }

            Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                << patches[patchi].name() << setprecision(3)
                << " " << setw(8)
                << returnReduce(patches[patchi].size(), sumOp<scalar>())
                << " " << setw(6) << layerParams.numLayers()[patchi]
                << " " << setw(8) << avgNearWallThickness
                << "  " << setw(8) << avgThickness
                << endl;
        }
        Info<< setprecision(oldPrecision) << endl;
    }
}


// Synchronize displacement among coupled patches.
void Foam::snappyLayerDriver::syncPatchDisplacement
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    //const fvMesh& mesh = meshRefiner.mesh();
    const labelList& meshPoints = pp.meshPoints();

    //label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax      // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        labelList syncPatchNLayers(patchNLayers);

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            minEqOp<label>(),
            labelMax            // null value
        );

        // Reset if differs
        // 1. take max
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Reset if differs
        // 2. take min
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }
        //nChangedTotal += nChanged;

        if (!returnReduceOr(nChanged))
        {
            break;
        }
    }

    //Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


// Calculate displacement vector for all patch points. Uses pointNormal.
// Checks that displaced patch point would be visible from all centres
// of the faces using it.
// extrudeStatus is both input and output and gives the status of each
// patch point.
void Foam::snappyLayerDriver::getPatchDisplacement
(
    const indirectPrimitivePatch& pp,
    const scalarField& thickness,
    const scalarField& minThickness,
    const scalarField& expansionRatio,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Determining displacement for added points"
        << " according to pointNormal ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();
    const pointField& localPoints = pp.localPoints();

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(PatchTools::pointNormals(mesh, pp));


    // Determine local length scale on patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off from same thickness everywhere (except where no extrusion)
    patchDisp = thickness*pointNormals;


    label nNoVisNormal = 0;
    label nExtrudeRemove = 0;


////XXXXXXXX
//    {
//        OBJstream twoStr
//        (
//            mesh.time().path()
//          / "twoFacePoints_"
//          + meshRefiner_.timeName()
//          + ".obj"
//        );
//        OBJstream multiStr
//        (
//            mesh.time().path()
//          / "multiFacePoints_"
//          + meshRefiner_.timeName()
//          + ".obj"
//        );
//        Pout<< "Writing points inbetween two faces on same cell to "
//            << twoStr.name() << endl;
//        Pout<< "Writing points inbetween three or more faces on same cell to "
//            << multiStr.name() << endl;
//        // Check whether inbetween patch faces on same cell
//        Map<labelList> cellToFaces;
//        forAll(pointNormals, patchPointi)
//        {
//            const labelList& pFaces = pointFaces[patchPointi];
//
//            cellToFaces.clear();
//            forAll(pFaces, pFacei)
//            {
//                const label patchFacei = pFaces[pFacei];
//                const label meshFacei = pp.addressing()[patchFacei];
//                const label celli = mesh.faceOwner()[meshFacei];
//
//                cellToFaces(celli).push_uniq(patchFacei);
//            }
//
//            forAllConstIters(cellToFaces, iter)
//            {
//                if (iter().size() == 2)
//                {
//                    twoStr.write(pp.localPoints()[patchPointi]);
//                }
//                else if (iter().size() > 2)
//                {
//                    multiStr.write(pp.localPoints()[patchPointi]);
//
//                    const scalar ratio =
//                        layerParameters::finalLayerThicknessRatio
//                        (
//                            patchNLayers[patchPointi],
//                            expansionRatio[patchPointi]
//                        );
//                    // Get thickness of cell next to bulk
//                    const vector finalDisp
//                    (
//                        ratio*patchDisp[patchPointi]
//                    );
//
//                    //Pout<< "** point:" << pp.localPoints()[patchPointi]
//                    //    << " on cell:" << iter.key()
//                    //    << " faces:" << iter()
//                    //    << " displacement was:" << patchDisp[patchPointi]
//                    //    << " ratio:" << ratio
//                    //    << " finalDispl:" << finalDisp;
//
//                    // Half this thickness
//                    patchDisp[patchPointi] -= 0.8*finalDisp;
//
//                    //Pout<< " new displacement:"
//                    //    << patchDisp[patchPointi] << endl;
//                }
//            }
//        }
//
//        Pout<< "Written " << multiStr.nVertices()
//            << " points inbetween three or more faces on same cell to "
//            << multiStr.name() << endl;
//    }
////XXXXXXXX


    // Check if no extrude possible.
    forAll(pointNormals, patchPointi)
    {
        label meshPointi = pp.meshPoints()[patchPointi];

        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Do not use unmarkExtrusion; forcibly set to zero extrusion.
            patchNLayers[patchPointi] = 0;
            patchDisp[patchPointi] = Zero;
        }
        else
        {
            // Get normal
            const vector& n = pointNormals[patchPointi];

            if (!meshTools::visNormal(n, faceNormals, pointFaces[patchPointi]))
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "No valid normal for point " << meshPointi
                        << ' ' << pp.points()[meshPointi]
                        << "; setting displacement to "
                        << patchDisp[patchPointi]
                        << endl;
                }

                extrudeStatus[patchPointi] = EXTRUDEREMOVE;
                nNoVisNormal++;
            }
        }
    }

    // At illegal points make displacement average of new neighbour positions
    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
        {
            point avg(Zero);
            label nPoints = 0;

            const labelList& pEdges = pp.pointEdges()[patchPointi];

            forAll(pEdges, i)
            {
                label edgei = pEdges[i];

                label otherPointi = pp.edges()[edgei].otherVertex(patchPointi);

                if (extrudeStatus[otherPointi] != NOEXTRUDE)
                {
                    avg += localPoints[otherPointi] + patchDisp[otherPointi];
                    nPoints++;
                }
            }

            if (nPoints > 0)
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "Displacement at illegal point "
                        << localPoints[patchPointi]
                        << " set to "
                        << (avg / nPoints - localPoints[patchPointi])
                        << endl;
                }

                patchDisp[patchPointi] =
                    avg / nPoints
                  - localPoints[patchPointi];

                nExtrudeRemove++;
            }
            else
            {
                // All surrounding points are not extruded. Leave patchDisp
                // intact.
            }
        }
    }

    Info<< "Detected " << returnReduce(nNoVisNormal, sumOp<label>())
        << " points with point normal pointing through faces." << nl
        << "Reset displacement at "
        << returnReduce(nExtrudeRemove, sumOp<label>())
        << " points to average of surrounding points." << endl;

    // Make sure displacement is equal on both sides of coupled patches.
    syncPatchDisplacement
    (
        mesh,
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    Info<< endl;
}


bool Foam::snappyLayerDriver::sameEdgeNeighbour
(
    const labelListList& globalEdgeFaces,
    const label myGlobalFacei,
    const label nbrGlobFacei,
    const label edgei
) const
{
    const labelList& eFaces = globalEdgeFaces[edgei];
    if (eFaces.size() == 2)
    {
        return edge(myGlobalFacei, nbrGlobFacei) == edge(eFaces[0], eFaces[1]);
    }

    return false;
}


void Foam::snappyLayerDriver::getVertexString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const label facei,
    const label edgei,
    const label myGlobFacei,
    const label nbrGlobFacei,
    DynamicList<label>& vertices
) const
{
    const labelList& fEdges = pp.faceEdges()[facei];
    label fp = fEdges.find(edgei);

    if (fp == -1)
    {
        FatalErrorInFunction
            << "problem." << abort(FatalError);
    }

    // Search back
    label startFp = fp;

    forAll(fEdges, i)
    {
        label prevFp = fEdges.rcIndex(startFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[prevFp]
            )
        )
        {
            break;
        }
        startFp = prevFp;
    }

    label endFp = fp;
    forAll(fEdges, i)
    {
        label nextFp = fEdges.fcIndex(endFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[nextFp]
            )
        )
        {
            break;
        }
        endFp = nextFp;
    }

    const face& f = pp.localFaces()[facei];
    vertices.clear();
    fp = startFp;
    while (fp != endFp)
    {
        vertices.append(f[fp]);
        fp = f.fcIndex(fp);
    }
    vertices.append(f[fp]);
    fp = f.fcIndex(fp);
    vertices.append(f[fp]);
}


// Truncates displacement
// - for all patchFaces in the faceset displacement gets set to zero
// - all displacement < minThickness gets set to zero
Foam::label Foam::snappyLayerDriver::truncateDisplacement
(
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    const faceSet& illegalPatchFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label nChanged = 0;

    const Map<label>& meshPointMap = pp.meshPointMap();

    for (const label facei : illegalPatchFaces)
    {
        if (mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Faceset " << illegalPatchFaces.name()
                << " contains internal face " << facei << nl
                << "It should only contain patch faces" << abort(FatalError);
        }

        const face& f = mesh.faces()[facei];


        forAll(f, fp)
        {
            const auto fnd = meshPointMap.cfind(f[fp]);
            if (fnd.good())
            {
                const label patchPointi = fnd.val();

                if (extrudeStatus[patchPointi] != NOEXTRUDE)
                {
                    unmarkExtrusion
                    (
                        patchPointi,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    nChanged++;
                }
            }
        }
    }

    forAll(patchDisp, patchPointi)
    {
        if (mag(patchDisp[patchPointi]) < minThickness[patchPointi])
        {
            if
            (
                unmarkExtrusion
                (
                    patchPointi,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nChanged++;
            }
        }
        else if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Make sure displacement is 0. Should already be so but ...
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }


    const faceList& localFaces = pp.localFaces();

    while (true)
    {
        syncPatchDisplacement
        (
            mesh,
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Pinch
        // ~~~~~

        // Make sure that a face doesn't have two non-consecutive areas
        // not extruded (e.g. quad where vertex 0 and 2 are not extruded
        // but 1 and 3 are) since this gives topological errors.

        label nPinched = 0;

        forAll(localFaces, i)
        {
            const face& localF = localFaces[i];

            // Count number of transitions from unsnapped to snapped.
            label nTrans = 0;

            extrudeMode prevMode = extrudeStatus[localF.prevLabel(0)];

            forAll(localF, fp)
            {
                extrudeMode fpMode = extrudeStatus[localF[fp]];

                if (prevMode == NOEXTRUDE && fpMode != NOEXTRUDE)
                {
                    nTrans++;
                }
                prevMode = fpMode;
            }

            if (nTrans > 1)
            {
                // Multiple pinches. Reset whole face as unextruded.
                if
                (
                    unmarkExtrusion
                    (
                        localF,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nPinched++;
                    nChanged++;
                }
            }
        }

        reduce(nPinched, sumOp<label>());

        Info<< "truncateDisplacement : Unextruded " << nPinched
            << " faces due to non-consecutive vertices being extruded." << endl;


        // Butterfly
        // ~~~~~~~~~

        // Make sure that a string of edges becomes a single face so
        // not a butterfly. Occasionally an 'edge' will have a single dangling
        // vertex due to face combining. These get extruded as a single face
        // (with a dangling vertex) so make sure this extrusion forms a single
        // shape.
        //  - continuous i.e. no butterfly:
        //      +     +
        //      |\   /|
        //      | \ / |
        //      +--+--+
        //  - extrudes from all but the endpoints i.e. no partial
        //    extrude
        //            +
        //           /|
        //          / |
        //      +--+--+
        // The common error topology is a pinch somewhere in the middle
        label nButterFly = 0;
        {
            DynamicList<label> stringedVerts;
            forAll(pp.edges(), edgei)
            {
                const labelList& globFaces = edgeGlobalFaces[edgei];

                if (globFaces.size() == 2)
                {
                    label myFacei = pp.edgeFaces()[edgei][0];
                    label myGlobalFacei = globalFaces.toGlobal
                    (
                        pp.addressing()[myFacei]
                    );
                    label nbrGlobalFacei =
                    (
                        globFaces[0] != myGlobalFacei
                      ? globFaces[0]
                      : globFaces[1]
                    );
                    getVertexString
                    (
                        pp,
                        edgeGlobalFaces,
                        myFacei,
                        edgei,
                        myGlobalFacei,
                        nbrGlobalFacei,
                        stringedVerts
                    );

                    if
                    (
                        extrudeStatus[stringedVerts[0]] != NOEXTRUDE
                     || extrudeStatus[stringedVerts.last()] != NOEXTRUDE
                    )
                    {
                        // Any pinch in the middle
                        bool pinch = false;
                        for (label i = 1; i < stringedVerts.size()-1; i++)
                        {
                            if (extrudeStatus[stringedVerts[i]] == NOEXTRUDE)
                            {
                                pinch = true;
                                break;
                            }
                        }
                        if (pinch)
                        {
                            forAll(stringedVerts, i)
                            {
                                if
                                (
                                    unmarkExtrusion
                                    (
                                        stringedVerts[i],
                                        patchDisp,
                                        patchNLayers,
                                        extrudeStatus
                                    )
                                )
                                {
                                    nButterFly++;
                                    nChanged++;
                                }
                            }
                        }
                    }
                }
            }
        }

        reduce(nButterFly, sumOp<label>());

        Info<< "truncateDisplacement : Unextruded " << nButterFly
            << " faces due to stringed edges with inconsistent extrusion."
            << endl;



        // Consistent number of layers
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Make sure that a face has consistent number of layers for all
        // its vertices.

        label nDiffering = 0;

        //forAll(localFaces, i)
        //{
        //    const face& localF = localFaces[i];
        //
        //    label numLayers = -1;
        //
        //    forAll(localF, fp)
        //    {
        //        if (patchNLayers[localF[fp]] > 0)
        //        {
        //            if (numLayers == -1)
        //            {
        //                numLayers = patchNLayers[localF[fp]];
        //            }
        //            else if (numLayers != patchNLayers[localF[fp]])
        //            {
        //                // Differing number of layers
        //                if
        //                (
        //                    unmarkExtrusion
        //                    (
        //                        localF,
        //                        patchDisp,
        //                        patchNLayers,
        //                        extrudeStatus
        //                    )
        //                )
        //                {
        //                    nDiffering++;
        //                    nChanged++;
        //                }
        //                break;
        //            }
        //        }
        //    }
        //}
        //
        //reduce(nDiffering, sumOp<label>());
        //
        //Info<< "truncateDisplacement : Unextruded " << nDiffering
        //    << " faces due to having differing number of layers." << endl;

        if (nPinched+nButterFly+nDiffering == 0)
        {
            break;
        }
    }

    return nChanged;
}


// Setup layer information (at points and faces) to modify mesh topology in
// regions where layer mesh terminates.
void Foam::snappyLayerDriver::setupLayerInfoTruncation
(
    const indirectPrimitivePatch& pp,
    const labelList& patchNLayers,
    const List<extrudeMode>& extrudeStatus,
    const label nBufferCellsNoExtrude,
    labelList& nPatchPointLayers,
    labelList& nPatchFaceLayers
) const
{
    Info<< nl << "Setting up information for layer truncation ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (nBufferCellsNoExtrude < 0)
    {
        Info<< nl << "Performing no layer truncation."
            << " nBufferCellsNoExtrude set to less than 0  ..." << endl;

        // Face layers if any point gets extruded
        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            forAll(f, fp)
            {
                const label nPointLayers = patchNLayers[f[fp]];
                if (nPointLayers > 0)
                {
                    if (nPatchFaceLayers[patchFacei] == -1)
                    {
                        nPatchFaceLayers[patchFacei] = nPointLayers;
                    }
                    else
                    {
                        nPatchFaceLayers[patchFacei] = min
                        (
                            nPatchFaceLayers[patchFacei],
                            nPointLayers
                        );
                    }
                }
            }
        }
        nPatchPointLayers = patchNLayers;

        // Set any unset patch face layers
        forAll(nPatchFaceLayers, patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = 0;
            }
        }
    }
    else
    {
        // Determine max point layers per face.
        labelList maxLevel(pp.size(), Zero);

        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            // find patch faces where layer terminates (i.e contains extrude
            // and noextrude points).

            bool noExtrude = false;
            label mLevel = 0;

            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    noExtrude = true;
                }
                mLevel = max(mLevel, patchNLayers[f[fp]]);
            }

            if (mLevel > 0)
            {
                // So one of the points is extruded. Check if all are extruded
                // or is a mix.

                if (noExtrude)
                {
                    nPatchFaceLayers[patchFacei] = 1;
                    maxLevel[patchFacei] = mLevel;
                }
                else
                {
                    maxLevel[patchFacei] = mLevel;
                }
            }
        }

        // We have the seed faces (faces with nPatchFaceLayers != maxLevel)
        // Now do a meshwave across the patch where we pick up neighbours
        // of seed faces.
        // Note: quite inefficient. Could probably be coded better.

        const labelListList& pointFaces = pp.pointFaces();

        label nLevels = gMax(patchNLayers);

        // flag neighbouring patch faces with number of layers to grow
        for (label ilevel = 1; ilevel < nLevels; ilevel++)
        {
            label nBuffer;

            if (ilevel == 1)
            {
                nBuffer = nBufferCellsNoExtrude - 1;
            }
            else
            {
                nBuffer = nBufferCellsNoExtrude;
            }

            for (label ibuffer = 0; ibuffer < nBuffer + 1; ibuffer++)
            {
                labelList tempCounter(nPatchFaceLayers);

                boolList foundNeighbour(pp.nPoints(), false);

                forAll(pp.meshPoints(), patchPointi)
                {
                    forAll(pointFaces[patchPointi], pointFacei)
                    {
                        label facei = pointFaces[patchPointi][pointFacei];

                        if
                        (
                            nPatchFaceLayers[facei] != -1
                         && maxLevel[facei] > 0
                        )
                        {
                            foundNeighbour[patchPointi] = true;
                            break;
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    pp.meshPoints(),
                    foundNeighbour,
                    orEqOp<bool>(),
                    false               // null value
                );

                forAll(pp.meshPoints(), patchPointi)
                {
                    if (foundNeighbour[patchPointi])
                    {
                        forAll(pointFaces[patchPointi], pointFacei)
                        {
                            label facei = pointFaces[patchPointi][pointFacei];
                            if
                            (
                                nPatchFaceLayers[facei] == -1
                             && maxLevel[facei] > 0
                             && ilevel < maxLevel[facei]
                            )
                            {
                                tempCounter[facei] = ilevel;
                            }
                        }
                    }
                }
                nPatchFaceLayers = tempCounter;
            }
        }

        forAll(pp.localFaces(), patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = maxLevel[patchFacei];
            }
        }

        forAll(pp.meshPoints(), patchPointi)
        {
            if (extrudeStatus[patchPointi] != NOEXTRUDE)
            {
                forAll(pointFaces[patchPointi], pointFacei)
                {
                    label face = pointFaces[patchPointi][pointFacei];
                    nPatchPointLayers[patchPointi] = max
                    (
                        nPatchPointLayers[patchPointi],
                        nPatchFaceLayers[face]
                    );
                }
            }
            else
            {
                nPatchPointLayers[patchPointi] = 0;
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPatchPointLayers,
            maxEqOp<label>(),
            label(0)        // null value
        );
    }
}


// Does any of the cells use a face from faces?
bool Foam::snappyLayerDriver::cellsUseFace
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelHashSet& faces
)
{
    forAll(cellLabels, i)
    {
        const cell& cFaces = mesh.cells()[cellLabels[i]];

        forAll(cFaces, cFacei)
        {
            if (faces.found(cFaces[cFacei]))
            {
                return true;
            }
        }
    }

    return false;
}


// Checks the newly added cells and locally unmarks points so they
// will not get extruded next time round. Returns global number of unmarked
// points (0 if all was fine)
Foam::label Foam::snappyLayerDriver::checkAndUnmark
(
    const addPatchCellLayer& addLayer,
    const dictionary& meshQualityDict,
    const bool additionalReporting,
    const List<labelPair>& baffles,
    const indirectPrimitivePatch& pp,
    const fvMesh& newMesh,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    // Check the resulting mesh for errors
    Info<< nl << "Checking mesh with layer ..." << endl;
    faceSet wrongFaces(newMesh, "wrongFaces", newMesh.nFaces()/1000);
    motionSmoother::checkMesh
    (
        false,
        newMesh,
        meshQualityDict,
        identity(newMesh.nFaces()),
        baffles,
        wrongFaces,
        false           // dryRun_
    );
    Info<< "Detected " << returnReduce(wrongFaces.size(), sumOp<label>())
        << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;

    // Undo local extrusion if
    // - any of the added cells in error

    label nChanged = 0;

    // Get all cells in the layer.
    labelListList addedCells
    (
        addPatchCellLayer::addedCells
        (
            newMesh,
            addLayer.layerFaces()
        )
    );

    // Check if any of the faces in error uses any face of an added cell
    // - if additionalReporting print the few remaining areas for ease of
    //   finding out where the problems are.

    const label nReportMax = 10;
    DynamicField<point> disabledFaceCentres(nReportMax);

    forAll(addedCells, oldPatchFacei)
    {
        // Get the cells (in newMesh labels) per old patch face (in mesh
        // labels)
        const labelList& fCells = addedCells[oldPatchFacei];

        if (cellsUseFace(newMesh, fCells, wrongFaces))
        {
            // Unmark points on old mesh
            if
            (
                unmarkExtrusion
                (
                    pp.localFaces()[oldPatchFacei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                if (additionalReporting && (nChanged < nReportMax))
                {
                    disabledFaceCentres.append
                    (
                        pp.faceCentres()[oldPatchFacei]
                    );
                }

                nChanged++;
            }
        }
    }


    label nChangedTotal = returnReduce(nChanged, sumOp<label>());

    if (additionalReporting)
    {
        // Limit the number of points to be printed so that
        // not too many points are reported when running in parallel
        // Not accurate, i.e. not always nReportMax points are written,
        // but this estimation avoid some communication here.
        // The important thing, however, is that when only a few faces
        // are disabled, their coordinates are printed, and this should be
        // the case
        label nReportLocal = nChanged;
        if (nChangedTotal > nReportMax)
        {
            nReportLocal = min
            (
                max(nChangedTotal / Pstream::nProcs(), 1),
                min
                (
                    nChanged,
                    max(nReportMax / Pstream::nProcs(), 1)
                )
            );
        }

        if (nReportLocal)
        {
            Pout<< "Checked mesh with layers. Disabled extrusion at " << endl;
            for (label i=0; i < nReportLocal; i++)
            {
                Pout<< "    " << disabledFaceCentres[i] << endl;
            }
        }

        label nReportTotal = returnReduce(nReportLocal, sumOp<label>());

        if (nReportTotal < nChangedTotal)
        {
            Info<< "Suppressed disabled extrusion message for other "
                << nChangedTotal - nReportTotal << " faces." << endl;
        }
    }

    return nChangedTotal;
}


//- Count global number of extruded faces
Foam::label Foam::snappyLayerDriver::countExtrusion
(
    const indirectPrimitivePatch& pp,
    const List<extrudeMode>& extrudeStatus
)
{
    // Count number of extruded patch faces
    label nExtruded = 0;
    {
        const faceList& localFaces = pp.localFaces();

        forAll(localFaces, i)
        {
            const face& localFace = localFaces[i];

            forAll(localFace, fp)
            {
                if (extrudeStatus[localFace[fp]] != NOEXTRUDE)
                {
                    nExtruded++;
                    break;
                }
            }
        }
    }

    return returnReduce(nExtruded, sumOp<label>());
}


Foam::List<Foam::labelPair> Foam::snappyLayerDriver::getBafflesOnAddedMesh
(
    const polyMesh& mesh,
    const labelList& newToOldFaces,
    const List<labelPair>& baffles
)
{
    // The problem is that the baffle faces are now inside the
    // mesh (addPatchCellLayer modifies original boundary faces and
    // adds new ones. So 2 pass:
    // - find the boundary face for all faces originating from baffle
    // - use the boundary face for the new baffles

    Map<label> baffleSet(4*baffles.size());
    forAll(baffles, bafflei)
    {
        baffleSet.insert(baffles[bafflei][0], bafflei);
        baffleSet.insert(baffles[bafflei][1], bafflei);
    }


    List<labelPair> newBaffles(baffles.size(), labelPair(-1, -1));
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label oldFacei = newToOldFaces[facei];

        const auto faceFnd = baffleSet.find(oldFacei);
        if (faceFnd.good())
        {
            label bafflei = faceFnd();
            labelPair& p = newBaffles[bafflei];
            if (p[0] == -1)
            {
                p[0] = facei;
            }
            else if (p[1] == -1)
            {
                p[1] = facei;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem:" << facei << " at:"
                    << mesh.faceCentres()[facei]
                    << " is on same baffle as " << p[0]
                    << " at:" << mesh.faceCentres()[p[0]]
                    << " and " << p[1]
                    << " at:" << mesh.faceCentres()[p[1]]
                    << exit(FatalError);
            }
        }
    }
    return newBaffles;
}


// Collect layer faces and layer cells into mesh fields for ease of handling
void Foam::snappyLayerDriver::getLayerCellsFaces
(
    const polyMesh& mesh,
    const addPatchCellLayer& addLayer,
    const scalarField& oldRealThickness,

    labelList& cellNLayers,
    scalarField& faceRealThickness
)
{
    cellNLayers.setSize(mesh.nCells());
    cellNLayers = 0;
    faceRealThickness.setSize(mesh.nFaces());
    faceRealThickness = 0;

    // Mark all faces in the layer
    const labelListList& layerFaces = addLayer.layerFaces();

    // Mark all cells in the layer.
    labelListList addedCells(addPatchCellLayer::addedCells(mesh, layerFaces));

    forAll(addedCells, oldPatchFacei)
    {
        const labelList& added = addedCells[oldPatchFacei];

        const labelList& layer = layerFaces[oldPatchFacei];

        if (layer.size())
        {
            // Leave out original internal face
            forAll(added, i)
            {
                cellNLayers[added[i]] = layer.size()-1;
            }
        }
    }

    forAll(layerFaces, oldPatchFacei)
    {
        const labelList& layer = layerFaces[oldPatchFacei];
        const scalar realThickness = oldRealThickness[oldPatchFacei];

        if (layer.size())
        {
            // Layer contains both original boundary face and new boundary
            // face so is nLayers+1. Leave out old internal face.
            for (label i = 1; i < layer.size(); i++)
            {
                faceRealThickness[layer[i]] = realThickness;
            }
        }
    }
}


void Foam::snappyLayerDriver::printLayerData
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& cellNLayers,
    const scalarField& faceWantedThickness,
    const scalarField& faceRealThickness,
    const layerParameters& layerParams
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const int oldPrecision = Info.stream().precision();

    // Find maximum length of a patch name, for a nicer output
    label maxPatchNameLen = 0;
    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        word patchName = pbm[patchi].name();
        maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
    }

    // Print some global mesh statistics
    meshRefiner_.printMeshInfo(false, "Mesh with layers", false);

    Info<< nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
        << setw(0) << " faces        layers        overall thickness" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << " "
        << setw(0) << "          target   mesh     [m]       [%]" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
        << setw(0) << " -----    -----    ----     ---       ---" << endl;


    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        const polyPatch& pp = pbm[patchi];

        label sumSize = pp.size();

        // Number of layers
        const labelList& faceCells = pp.faceCells();
        label sumNLayers = 0;
        forAll(faceCells, i)
        {
            sumNLayers += cellNLayers[faceCells[i]];
        }

        // Thickness
        scalarField::subField patchWanted = pbm[patchi].patchSlice
        (
            faceWantedThickness
        );
        scalarField::subField patchReal = pbm[patchi].patchSlice
        (
            faceRealThickness
        );

        scalar sumRealThickness = sum(patchReal);
        scalar sumFraction = 0;
        forAll(patchReal, i)
        {
            if (patchWanted[i] > VSMALL)
            {
                sumFraction += (patchReal[i]/patchWanted[i]);
            }
        }


        reduce(sumSize, sumOp<label>());
        reduce(sumNLayers, sumOp<label>());
        reduce(sumRealThickness, sumOp<scalar>());
        reduce(sumFraction, sumOp<scalar>());


        scalar avgLayers = 0;
        scalar avgReal = 0;
        scalar avgFraction = 0;
        if (sumSize > 0)
        {
            avgLayers = scalar(sumNLayers)/sumSize;
            avgReal = sumRealThickness/sumSize;
            avgFraction = sumFraction/sumSize;
        }

        Info<< setf(ios_base::left) << setw(maxPatchNameLen)
            << pbm[patchi].name() << setprecision(3)
            << " " << setw(8) << sumSize
            << " " << setw(8) << layerParams.numLayers()[patchi]
            << " " << setw(8) << avgLayers
            << " " << setw(8) << avgReal
            << "  " << setw(8) << 100*avgFraction
            << endl;
    }
    Info<< setprecision(oldPrecision) << endl;
}


bool Foam::snappyLayerDriver::writeLayerSets
(
    const fvMesh& mesh,
    const labelList& cellNLayers,
    const scalarField& faceRealThickness
) const
{
    bool allOk = true;
    {
        label nAdded = 0;
        forAll(cellNLayers, celli)
        {
            if (cellNLayers[celli] > 0)
            {
                nAdded++;
            }
        }
        cellSet addedCellSet(mesh, "addedCells", nAdded);
        forAll(cellNLayers, celli)
        {
            if (cellNLayers[celli] > 0)
            {
                addedCellSet.insert(celli);
            }
        }
        addedCellSet.instance() = meshRefiner_.timeName();
        Info<< "Writing "
            << returnReduce(addedCellSet.size(), sumOp<label>())
            << " added cells to cellSet "
            << addedCellSet.name() << endl;
        bool ok = addedCellSet.write();
        allOk = allOk && ok;
    }
    {
        label nAdded = 0;
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (faceRealThickness[facei] > 0)
            {
                nAdded++;
            }
        }

        faceSet layerFacesSet(mesh, "layerFaces", nAdded);
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (faceRealThickness[facei] > 0)
            {
                layerFacesSet.insert(facei);
            }
        }
        layerFacesSet.instance() = meshRefiner_.timeName();
        Info<< "Writing "
            << returnReduce(layerFacesSet.size(), sumOp<label>())
            << " faces inside added layer to faceSet "
            << layerFacesSet.name() << endl;
        bool ok = layerFacesSet.write();
        allOk = allOk && ok;
    }
    return allOk;
}


bool Foam::snappyLayerDriver::writeLayerData
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& cellNLayers,
    const scalarField& faceWantedThickness,
    const scalarField& faceRealThickness
) const
{
    bool allOk = true;

    if (meshRefinement::writeLevel() & meshRefinement::WRITELAYERSETS)
    {
        bool ok = writeLayerSets(mesh, cellNLayers, faceRealThickness);
        allOk = allOk && ok;
    }

    if (meshRefinement::writeLevel() & meshRefinement::WRITELAYERFIELDS)
    {
        Info<< nl << "Writing fields with layer information:" << incrIndent
            << endl;
        {
            volScalarField fld
            (
                IOobject
                (
                    "nSurfaceLayers",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimless, Zero),
                fixedValueFvPatchScalarField::typeName
            );
            forAll(fld, celli)
            {
                fld[celli] = cellNLayers[celli];
            }
            volScalarField::Boundary& fldBf = fld.boundaryFieldRef();

            const polyBoundaryMesh& pbm = mesh.boundaryMesh();
            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];
                const polyPatch& pp = pbm[patchi];
                const labelList& faceCells = pp.faceCells();
                scalarField pfld(faceCells.size());
                forAll(faceCells, i)
                {
                    pfld[i] = cellNLayers[faceCells[i]];
                }
                fldBf[patchi] == pfld;
            }
            Info<< indent << fld.name() << "    : actual number of layers"
                << endl;
            bool ok = fld.write();
            allOk = allOk && ok;
        }
        {
            volScalarField fld
            (
                IOobject
                (
                    "thickness",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimless, Zero),
                fixedValueFvPatchScalarField::typeName
            );
            volScalarField::Boundary& fldBf = fld.boundaryFieldRef();
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();
            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];
                fldBf[patchi] == pbm[patchi].patchSlice(faceRealThickness);
            }
            Info<< indent << fld.name() << "         : overall layer thickness"
                << endl;
            bool ok = fld.write();
            allOk = allOk && ok;
        }
        {
            volScalarField fld
            (
                IOobject
                (
                    "thicknessFraction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimless, Zero),
                fixedValueFvPatchScalarField::typeName
            );
            volScalarField::Boundary& fldBf = fld.boundaryFieldRef();
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();
            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];

                scalarField::subField patchWanted = pbm[patchi].patchSlice
                (
                    faceWantedThickness
                );
                scalarField::subField patchReal = pbm[patchi].patchSlice
                (
                    faceRealThickness
                );

                // Convert patchReal to relative thickness
                scalarField pfld(patchReal.size(), Zero);
                forAll(patchReal, i)
                {
                    if (patchWanted[i] > VSMALL)
                    {
                        pfld[i] = patchReal[i]/patchWanted[i];
                    }
                }

                fldBf[patchi] == pfld;
            }
            Info<< indent << fld.name()
                << " : overall layer thickness (fraction"
                << " of desired thickness)" << endl;
            bool ok = fld.write();
            allOk = allOk && ok;
        }
        Info<< decrIndent<< endl;
    }

    return allOk;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::snappyLayerDriver::dupFaceZonePoints
(
    meshRefinement& meshRefiner,
    const labelList& patchIDs,  // patch indices
    const labelList& numLayers, // number of layers per patch
    List<labelPair> baffles,    // pairs of baffles (input & updated)
    labelList& pointToMaster    // -1 or index of original point (duplicated
                                // point)
)
{
    fvMesh& mesh = meshRefiner.mesh();

    // Check outside of baffles for non-manifoldness

    // Points that are candidates for duplication
    labelList candidatePoints;
    {
        // Do full analysis to see if we need to extrude points
        // so have to duplicate them
        autoPtr<indirectPrimitivePatch> pp
        (
            meshRefinement::makePatch
            (
                mesh,
                patchIDs
            )
        );

        // Displacement for all pp.localPoints. Set to large value to
        // avoid truncation in syncPatchDisplacement because of
        // minThickness.
        vectorField patchDisp(pp().nPoints(), vector(GREAT, GREAT, GREAT));
        labelList patchNLayers(pp().nPoints(), Zero);
        label nIdealTotAddedCells = 0;
        List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);
        // Get number of layers per point from number of layers per patch
        setNumLayers
        (
            meshRefiner,
            numLayers,              // per patch the num layers
            patchIDs,               // patches that are being moved
            *pp,                    // indirectpatch for all faces moving

            patchNLayers,
            extrudeStatus,
            nIdealTotAddedCells
        );
        // Make sure displacement is equal on both sides of coupled patches.
        // Note that we explicitly disable the minThickness truncation
        // of the patchDisp here.
        syncPatchDisplacement
        (
            mesh,
            *pp,
            scalarField(patchDisp.size(), Zero), //minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Do duplication only if all patch points decide to extrude. Ignore
        // contribution from non-patch points. Note that we need to
        // apply this to all mesh points
        labelList minPatchState(mesh.nPoints(), labelMax);
        forAll(extrudeStatus, patchPointi)
        {
            label pointi = pp().meshPoints()[patchPointi];
            minPatchState[pointi] = extrudeStatus[patchPointi];
        }

        syncTools::syncPointList
        (
            mesh,
            minPatchState,
            minEqOp<label>(),   // combine op
            labelMax            // null value
        );

        // So now minPatchState:
        // - labelMax on non-patch points
        // - NOEXTRUDE if any patch point was not extruded
        // - EXTRUDE or EXTRUDEREMOVE if all patch points are extruded/
        //   extrudeRemove.

        label n = 0;
        forAll(minPatchState, pointi)
        {
            label state = minPatchState[pointi];
            if (state == EXTRUDE || state == EXTRUDEREMOVE)
            {
                n++;
            }
        }
        candidatePoints.setSize(n);
        n = 0;
        forAll(minPatchState, pointi)
        {
            label state = minPatchState[pointi];
            if (state == EXTRUDE || state == EXTRUDEREMOVE)
            {
                candidatePoints[n++] = pointi;
            }
        }
    }

    // Not duplicate points on either side of baffles that don't get any
    // layers
    labelPairList nonDupBaffles;

    {
        // faceZones that are not being duplicated
        DynamicList<label> nonDupZones(mesh.faceZones().size());

        labelHashSet layerIDs(patchIDs);
        forAll(mesh.faceZones(), zonei)
        {
            label mpi, spi;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner.getFaceZoneInfo
            (
                mesh.faceZones()[zonei].name(),
                mpi,
                spi,
                fzType
            );
            if (hasInfo && !layerIDs.found(mpi) && !layerIDs.found(spi))
            {
                nonDupZones.append(zonei);
            }
        }
        nonDupBaffles = meshRefinement::subsetBaffles
        (
            mesh,
            nonDupZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        );
    }


    const localPointRegion regionSide(mesh, nonDupBaffles, candidatePoints);

    autoPtr<mapPolyMesh> map = meshRefiner.dupNonManifoldPoints
    (
        regionSide
    );

    if (map)
    {
        // Store point duplication
        pointToMaster.setSize(mesh.nPoints(), -1);

        const labelList& pointMap = map().pointMap();
        const labelList& reversePointMap = map().reversePointMap();

        forAll(pointMap, pointi)
        {
            label oldPointi = pointMap[pointi];
            label newMasterPointi = reversePointMap[oldPointi];

            if (newMasterPointi != pointi)
            {
                // Found slave. Mark both master and slave
                pointToMaster[pointi] = newMasterPointi;
                pointToMaster[newMasterPointi] = newMasterPointi;
            }
        }

        // Update baffle numbering
        {
            const labelList& reverseFaceMap = map().reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }


        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing point-duplicate mesh to time "
                << meshRefiner.timeName() << endl;

            meshRefiner.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner.timeName()
            );

            OBJstream str
            (
                mesh.time().path()
              / "duplicatePoints_"
              + meshRefiner.timeName()
              + ".obj"
            );
            Info<< "Writing point-duplicates to " << str.name() << endl;
            const pointField& p = mesh.points();
            forAll(pointMap, pointi)
            {
                label newMasteri = reversePointMap[pointMap[pointi]];

                if (newMasteri != pointi)
                {
                    str.writeLine(p[pointi], p[newMasteri]);
                }
            }
        }
    }
    return map;
}


void Foam::snappyLayerDriver::mergeFaceZonePoints
(
    const labelList& pointToMaster, // -1 or index of duplicated point
    labelList& cellNLayers,
    scalarField& faceRealThickness,
    scalarField& faceWantedThickness
)
{
    // Near opposite of dupFaceZonePoints : merge points and baffles introduced
    // for internal faceZones

    fvMesh& mesh = meshRefiner_.mesh();

    // Count duplicate points
    label nPointPairs = 0;
    forAll(pointToMaster, pointi)
    {
        label otherPointi = pointToMaster[pointi];
        if (otherPointi != -1)
        {
            nPointPairs++;
        }
    }
    reduce(nPointPairs, sumOp<label>());
    if (nPointPairs > 0)
    {
        // Merge any duplicated points
        Info<< "Merging " << nPointPairs << " duplicated points ..." << endl;

        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            OBJstream str
            (
                mesh.time().path()
              / "mergePoints_"
              + meshRefiner_.timeName()
              + ".obj"
            );
            Info<< "Points to be merged to " << str.name() << endl;
            forAll(pointToMaster, pointi)
            {
                label otherPointi = pointToMaster[pointi];
                if (otherPointi != -1)
                {
                    const point& pt = mesh.points()[pointi];
                    const point& otherPt = mesh.points()[otherPointi];
                    str.writeLine(pt, otherPt);
                }
            }
        }


        autoPtr<mapPolyMesh> map = meshRefiner_.mergePoints(pointToMaster);
        if (map)
        {
            inplaceReorder(map().reverseCellMap(), cellNLayers);

            const labelList& reverseFaceMap = map().reverseFaceMap();
            inplaceReorder(reverseFaceMap, faceWantedThickness);
            inplaceReorder(reverseFaceMap, faceRealThickness);

            Info<< "Merged points in = "
                << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
        }
    }

    if (mesh.faceZones().size() > 0)
    {
        // Merge any baffles
        Info<< "Converting baffles back into zoned faces ..."
            << endl;

        autoPtr<mapPolyMesh> map = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false   // baffle zones
        );
        if (map)
        {
            inplaceReorder(map().reverseCellMap(), cellNLayers);

            const labelList& faceMap = map().faceMap();

            // Make sure to keep the max since on two patches only one has
            // layers.
            scalarField newFaceRealThickness(mesh.nFaces(), Zero);
            scalarField newFaceWantedThickness(mesh.nFaces(), Zero);
            forAll(newFaceRealThickness, facei)
            {
                label oldFacei = faceMap[facei];
                if (oldFacei >= 0)
                {
                    scalar& realThick = newFaceRealThickness[facei];
                    realThick = max(realThick, faceRealThickness[oldFacei]);
                    scalar& wanted = newFaceWantedThickness[facei];
                    wanted = max(wanted, faceWantedThickness[oldFacei]);
                }
            }
            faceRealThickness.transfer(newFaceRealThickness);
            faceWantedThickness.transfer(newFaceWantedThickness);
        }

        Info<< "Converted baffles in = "
            << meshRefiner_.mesh().time().cpuTimeIncrement()
            << " s\n" << nl << endl;
    }
}


Foam::label Foam::snappyLayerDriver::setPointNumLayers
(
    const layerParameters& layerParams,

    const labelList& numLayers,
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    const labelListList& edgeGlobalFaces,

    vectorField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    patchDisp.setSize(pp.nPoints());
    patchDisp = vector(GREAT, GREAT, GREAT);

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    patchNLayers.setSize(pp.nPoints());
    patchNLayers = Zero;

    // Ideal number of cells added
    label nIdealTotAddedCells = 0;

    // Whether to add edge for all pp.localPoints.
    extrudeStatus.setSize(pp.nPoints());
    extrudeStatus = EXTRUDE;


    // Get number of layers per point from number of layers per patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    setNumLayers
    (
        meshRefiner_,
        numLayers,                  // per patch the num layers
        patchIDs,                   // patches that are being moved
        pp,                         // indirectpatch for all faces moving

        patchNLayers,
        extrudeStatus,
        nIdealTotAddedCells
    );

    // Precalculate mesh edge labels for patch edges
    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));


    // Disable extrusion on split strings of common points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    handleNonStringConnected
    (
        pp,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );


    // Disable extrusion on non-manifold points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    handleNonManifolds
    (
        pp,
        meshEdges,
        edgeGlobalFaces,

        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    // Disable extrusion on feature angles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    handleFeatureAngle
    (
        pp,
        meshEdges,
        layerParams.featureAngle(),

        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    // Disable extrusion on warped faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // It is hard to calculate some length scale if not in relative
    // mode so disable this check.
    if (!layerParams.relativeSizes().found(false))
    {
        // Undistorted edge length
        const scalar edge0Len =
            meshRefiner_.meshCutter().level0EdgeLength();
        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

        handleWarpedFaces
        (
            pp,
            layerParams.maxFaceThicknessRatio(),
            layerParams.relativeSizes(),
            edge0Len,
            cellLevel,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );
    }

    //// Disable extrusion on cells with multiple patch faces
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //handleMultiplePatchFaces
    //(
    //    pp,
    //
    //    patchDisp,
    //    patchNLayers,
    //    extrudeStatus
    //);

    addProfiling(grow, "snappyHexMesh::layers::grow");

    // Grow out region of non-extrusion
    for (label i = 0; i < layerParams.nGrow(); i++)
    {
        growNoExtrusion
        (
            pp,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );
    }
    return nIdealTotAddedCells;
}


Foam::autoPtr<Foam::externalDisplacementMeshMover>
Foam::snappyLayerDriver::makeMeshMover
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& internalFaceZones,
    const scalarIOField& minThickness,
    pointVectorField& displacement
) const
{
    // Allocate run-time selectable mesh mover

    fvMesh& mesh = meshRefiner_.mesh();


    // Set up controls for meshMover
    dictionary combinedDict(layerParams.dict());
    // Add mesh quality constraints
    combinedDict.merge(motionDict);
    // Where to get minThickness from
    combinedDict.add("minThicknessName", minThickness.name());

    const List<labelPair> internalBaffles
    (
        meshRefinement::subsetBaffles
        (
            mesh,
            internalFaceZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        )
    );

    // Take over patchDisp as boundary conditions on displacement
    // pointVectorField
    autoPtr<Foam::externalDisplacementMeshMover> medialAxisMoverPtr
    (
        externalDisplacementMeshMover::New
        (
            layerParams.meshShrinker(),
            combinedDict,
            internalBaffles,
            displacement
        )
    );


    if (dryRun_)
    {
        string errorMsg(FatalError.message());
        string IOerrorMsg(FatalIOError.message());

        if (errorMsg.size() || IOerrorMsg.size())
        {
            //errorMsg = "[dryRun] " + errorMsg;
            //errorMsg.replaceAll("\n", "\n[dryRun] ");
            //IOerrorMsg = "[dryRun] " + IOerrorMsg;
            //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

            IOWarningInFunction(combinedDict)
                << nl
                << "Missing/incorrect required dictionary entries:"
                << nl << nl
                << IOerrorMsg.c_str() << nl
                << errorMsg.c_str() << nl << nl
                << "Exiting dry-run" << nl << endl;

            if (UPstream::parRun())
            {
                Perr<< "\nFOAM parallel run exiting\n" << endl;
                UPstream::exit(0);
            }
            else
            {
                Perr<< "\nFOAM exiting\n" << endl;
                std::exit(0);
            }
        }
    }
    return medialAxisMoverPtr;
}


void Foam::snappyLayerDriver::addLayers
(
    const layerParameters& layerParams,
    const label nLayerIter,

    // Mesh quality provision
    const dictionary& motionDict,
    const label nRelaxedIter,
    const label nAllowableErrors,

    const labelList& patchIDs,
    const labelList& internalFaceZones,
    const List<labelPair>& baffles,
    const labelList& numLayers,
    const label nIdealTotAddedCells,

    const globalIndex& globalFaces,
    indirectPrimitivePatch& pp,

    const labelListList& edgeGlobalFaces,
    const labelList& edgePatchID,
    const labelList& edgeZoneID,
    const boolList& edgeFlip,
    const labelList& inflateFaceID,

    const scalarField& thickness,
    const scalarIOField& minThickness,
    const scalarField& expansionRatio,

    // Displacement for all pp.localPoints. Set to large value so does
    // not trigger the minThickness truncation (see syncPatchDisplacement below)
    vectorField& patchDisp,

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList& patchNLayers,

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode>& extrudeStatus,


    polyTopoChange& savedMeshMod,


    // Per cell 0 or number of layers in the cell column it is part of
    labelList& cellNLayers,
    // Per face actual overall layer thickness
    scalarField& faceRealThickness
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Overall displacement field
    pointVectorField displacement
    (
        makeLayerDisplacementField
        (
            pointMesh::New(mesh),
            numLayers
        )
    );

    // Allocate run-time selectable mesh mover
    autoPtr<externalDisplacementMeshMover> medialAxisMoverPtr
    (
        makeMeshMover
        (
            layerParams,
            motionDict,
            internalFaceZones,
            minThickness,
            displacement
        )
    );


    // Saved old points
    const pointField oldPoints(mesh.points());

    for (label iteration = 0; iteration < nLayerIter; iteration++)
    {
        Info<< nl
            << "Layer addition iteration " << iteration << nl
            << "--------------------------" << endl;


        // Unset the extrusion at the pp.
        const dictionary& meshQualityDict =
        (
            iteration < nRelaxedIter
          ? motionDict
          : motionDict.subDict("relaxed")
        );

        if (iteration >= nRelaxedIter)
        {
            Info<< "Switched to relaxed meshQuality constraints." << endl;
        }



        // Make sure displacement is equal on both sides of coupled patches.
        // Note that this also does the patchDisp < minThickness truncation
        // so for the first pass make sure the patchDisp is larger than
        // that.
        syncPatchDisplacement
        (
            mesh,
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Displacement acc. to pointnormals
        getPatchDisplacement
        (
            pp,
            thickness,
            minThickness,
            expansionRatio,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Shrink mesh by displacement value first.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            const pointField oldPatchPos(pp.localPoints());

            // We have patchDisp which is the outwards pointing
            // extrusion distance. Convert into an inwards pointing
            // shrink distance
            patchDisp = -patchDisp;

            // Take over patchDisp into pointDisplacement field and
            // adjust both for multi-patch constraints
            motionSmootherAlgo::setDisplacement
            (
                patchIDs,
                pp,
                patchDisp,
                displacement
            );


            // Move mesh
            // ~~~~~~~~~

            // Set up controls for meshMover
            dictionary combinedDict(layerParams.dict());
            // Add standard quality constraints
            combinedDict.merge(motionDict);
            // Add relaxed constraints (overrides standard ones)
            combinedDict.merge(meshQualityDict);
            // Where to get minThickness from
            combinedDict.add("minThicknessName", minThickness.name());

            labelList checkFaces(identity(mesh.nFaces()));
            medialAxisMoverPtr().move
            (
                combinedDict,
                nAllowableErrors,
                checkFaces
            );

            pp.movePoints(mesh.points());

            // Unset any moving state
            mesh.moving(false);

            // Update patchDisp (since not all might have been honoured)
            patchDisp = oldPatchPos - pp.localPoints();
        }

        // Truncate displacements that are too small (this will do internal
        // ones, coupled ones have already been truncated by
        // syncPatchDisplacement)
        faceSet dummySet(mesh, "wrongPatchFaces", 0);
        truncateDisplacement
        (
            globalFaces,
            edgeGlobalFaces,
            pp,
            minThickness,
            dummySet,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Dump to .obj file for debugging.
        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            dumpDisplacement
            (
                mesh.time().path()/"layer_" + meshRefiner_.timeName(),
                pp,
                patchDisp,
                extrudeStatus
            );

            const_cast<Time&>(mesh.time())++;
            Info<< "Writing shrunk mesh to time "
                << meshRefiner_.timeName() << endl;

            // See comment in snappySnapDriver why we should not remove
            // meshPhi using mesh.clearOut().

            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }


        // Mesh topo change engine. Insert current mesh.
        polyTopoChange meshMod(mesh);

        // Grow layer of cells on to patch. Handles zero sized displacement.
        addPatchCellLayer addLayer(mesh);

        // Determine per point/per face number of layers to extrude. Also
        // handles the slow termination of layers when going switching
        // layers

        labelList nPatchPointLayers(pp.nPoints(), -1);
        labelList nPatchFaceLayers(pp.size(), -1);
        setupLayerInfoTruncation
        (
            pp,
            patchNLayers,
            extrudeStatus,
            layerParams.nBufferCellsNoExtrude(),
            nPatchPointLayers,
            nPatchFaceLayers
        );

        // Calculate displacement for final layer for addPatchLayer.
        // (layer of cells next to the original mesh)
        vectorField finalDisp(patchNLayers.size(), Zero);

        forAll(nPatchPointLayers, i)
        {
            scalar ratio = layerParameters::finalLayerThicknessRatio
            (
                nPatchPointLayers[i],
                expansionRatio[i]
            );
            finalDisp[i] = ratio*patchDisp[i];
        }


        const scalarField invExpansionRatio(1.0/expansionRatio);

        // Add topo regardless of whether extrudeStatus is extruderemove.
        // Not add layer if patchDisp is zero.
        addLayer.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            invExpansionRatio,
            pp,
            bitSet(pp.size()),    // no flip

            edgePatchID,    // boundary patch for extruded boundary edges
            edgeZoneID,     // zone for extruded edges
            edgeFlip,
            inflateFaceID,


            labelList(0),   // exposed patchIDs, not used for adding layers
            nPatchFaceLayers,   // layers per face
            nPatchPointLayers,  // layers per point
            finalDisp,      // thickness of layer nearest internal mesh
            meshMod
        );

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Compact storage
        meshMod.shrink();

        // Store mesh changes for if mesh is correct.
        savedMeshMod = meshMod;


        // With the stored topo changes we create a new mesh so we can
        // undo if necessary.

        autoPtr<fvMesh> newMeshPtr;
        autoPtr<mapPolyMesh> mapPtr = meshMod.makeMesh
        (
            newMeshPtr,
            IOobject
            (
                //mesh.name()+"_layer",
                mesh.name(),
                static_cast<polyMesh&>(mesh).instance(),
                mesh.time(),  // register with runTime
                IOobject::READ_IF_PRESENT,  // read fv* if present
                static_cast<polyMesh&>(mesh).writeOpt()
            ),              // io params from original mesh but new name
            mesh,           // original mesh
            true            // parallel sync
        );
        fvMesh& newMesh = *newMeshPtr;
        mapPolyMesh& map = *mapPtr;

        // Get timing, but more importantly get memory information
        addProfiling(grow, "snappyHexMesh::layers::updateMesh");

        //?necessary? Update fields
        newMesh.updateMesh(map);

        // Move mesh if in inflation mode
        if (map.hasMotionPoints())
        {
            newMesh.movePoints(map.preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            newMesh.clearOut();
        }

        newMesh.setInstance(meshRefiner_.timeName());

        // Update numbering on addLayer:
        // - cell/point labels to be newMesh.
        // - patchFaces to remain in oldMesh order.
        addLayer.updateMesh
        (
            map,
            identity(pp.size()),
            identity(pp.nPoints())
        );

        // Collect layer faces and cells for outside loop.
        getLayerCellsFaces
        (
            newMesh,
            addLayer,
            avgPointData(pp, mag(patchDisp))(), // current thickness

            cellNLayers,
            faceRealThickness
        );


        // Count number of added cells
        label nAddedCells = 0;
        forAll(cellNLayers, celli)
        {
            if (cellNLayers[celli] > 0)
            {
                nAddedCells++;
            }
        }


        if (debug&meshRefinement::MESH)
        {
            Info<< "Writing layer mesh to time " << meshRefiner_.timeName()
                << endl;
            newMesh.write();
            writeLayerSets(newMesh, cellNLayers, faceRealThickness);

            // Reset the instance of the original mesh so next iteration
            // it dumps a complete mesh. This is just so that the inbetween
            // newMesh does not upset e.g. paraFoam cycling through the
            // times.
            mesh.setInstance(meshRefiner_.timeName());
        }


        //- Get baffles in newMesh numbering. Note that we cannot detect
        //  baffles here since the points are duplicated
        List<labelPair> internalBaffles;
        {
            // From old mesh face to corresponding newMesh boundary face
            labelList meshToNewMesh(mesh.nFaces(), -1);
            for
            (
                label facei = newMesh.nInternalFaces();
                facei < newMesh.nFaces();
                facei++
            )
            {
                label newMeshFacei = map.faceMap()[facei];
                if (newMeshFacei != -1)
                {
                    meshToNewMesh[newMeshFacei] = facei;
                }
            }

            List<labelPair> newMeshBaffles(baffles.size());
            label newi = 0;
            forAll(baffles, i)
            {
                const labelPair& p = baffles[i];
                labelPair newMeshBaffle
                (
                    meshToNewMesh[p[0]],
                    meshToNewMesh[p[1]]
                );
                if (newMeshBaffle[0] != -1 && newMeshBaffle[1] != -1)
                {
                    newMeshBaffles[newi++] = newMeshBaffle;
                }
            }
            newMeshBaffles.setSize(newi);

            internalBaffles = meshRefinement::subsetBaffles
            (
                newMesh,
                internalFaceZones,
                newMeshBaffles
            );

            Info<< "Detected "
                << returnReduce(internalBaffles.size(), sumOp<label>())
                << " baffles across faceZones of type internal" << nl
                << endl;
        }

        label nTotChanged = checkAndUnmark
        (
            addLayer,
            meshQualityDict,
            layerParams.additionalReporting(),
            internalBaffles,
            pp,
            newMesh,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        label nTotExtruded = countExtrusion(pp, extrudeStatus);
        label nTotFaces = returnReduce(pp.size(), sumOp<label>());
        label nTotAddedCells = returnReduce(nAddedCells, sumOp<label>());

        Info<< "Extruding " << nTotExtruded
            << " out of " << nTotFaces
            << " faces (" << 100.0*nTotExtruded/nTotFaces << "%)."
            << " Removed extrusion at " << nTotChanged << " faces."
            << endl
            << "Added " << nTotAddedCells << " out of "
            << nIdealTotAddedCells
            << " cells (" << 100.0*nTotAddedCells/nIdealTotAddedCells
            << "%)." << endl;

        if (nTotChanged == 0)
        {
            break;
        }

        // Reset mesh points and start again
        mesh.movePoints(oldPoints);
        pp.movePoints(mesh.points());
        medialAxisMoverPtr().movePoints(mesh.points());

        // Unset any moving state
        mesh.moving(false);

        // Grow out region of non-extrusion
        for (label i = 0; i < layerParams.nGrow(); i++)
        {
            growNoExtrusion
            (
                pp,
                patchDisp,
                patchNLayers,
                extrudeStatus
            );
        }

        Info<< endl;
    }
}


void Foam::snappyLayerDriver::mapFaceZonePoints
(
    meshRefinement& meshRefiner,
    const mapPolyMesh& map,
    labelPairList& baffles,
    labelList& pointToMaster
)
{
    fvMesh& mesh = meshRefiner.mesh();

    // Use geometric detection of points-to-be-merged
    //  - detect any boundary face created from a duplicated face (=baffle)
    //  - on these mark any point created from a duplicated point
    if (returnReduceOr(pointToMaster.size()))
    {
        // Estimate number of points-to-be-merged
        DynamicList<label> candidates(baffles.size()*4);

        // The problem is that all the internal layer faces also
        // have reverseFaceMap pointing to the old baffle face. So instead
        // loop over all the boundary faces and see which pair of new boundary
        // faces corresponds to the old baffles.


        // Mark whether old face was on baffle
        Map<label> oldFaceToBaffle(2*baffles.size());
        forAll(baffles, i)
        {
            const labelPair& baffle = baffles[i];
            oldFaceToBaffle.insert(baffle[0], i);
            oldFaceToBaffle.insert(baffle[1], i);
        }

        labelPairList newBaffles(baffles.size(), labelPair(-1, -1));

        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            const label oldFacei = map.faceMap()[facei];
            const auto iter = oldFaceToBaffle.find(oldFacei);
            if (oldFacei != -1 && iter.good())
            {
                const label bafflei = iter();
                auto& newBaffle = newBaffles[bafflei];
                if (newBaffle[0] == -1)
                {
                    newBaffle[0] = facei;
                }
                else if (newBaffle[1] == -1)
                {
                    newBaffle[1] = facei;
                }
                else
                {
                    FatalErrorInFunction << "face:" << facei
                        << " at:" << mesh.faceCentres()[facei]
                        << " already maps to baffle faces:"
                        << newBaffle[0]
                        << " at:" << mesh.faceCentres()[newBaffle[0]]
                        << " and " << newBaffle[1]
                        << " at:" << mesh.faceCentres()[newBaffle[1]]
                        << exit(FatalError);
                }

                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label oldPointi = map.pointMap()[pointi];

                    if (pointToMaster[oldPointi] != -1)
                    {
                        candidates.append(pointi);
                    }
                }
            }
        }


        // Compact newBaffles
        {
            label n = 0;
            forAll(newBaffles, i)
            {
                const labelPair& newBaffle = newBaffles[i];
                if (newBaffle[0] != -1 && newBaffle[1] != -1)
                {
                    newBaffles[n++] = newBaffle;
                }
            }

            newBaffles.setSize(n);
            baffles.transfer(newBaffles);
        }


        // Do geometric merge. Ideally we'd like to use a topological
        // merge but we've thrown away all layer-wise addressing when
        // throwing away addPatchCellLayer engine. Also the addressing
        // is extremely complicated. There is no problem with merging
        // too many points; the problem would be if merging baffles.
        // Trust mergeZoneBaffles to do sufficient checks.
        labelList oldToNew;
        label nNew = Foam::mergePoints
        (
            UIndirectList<point>(mesh.points(), candidates),
            meshRefiner.mergeDistance(),
            false,
            oldToNew
        );

        // Extract points to be merged (i.e. multiple points originating
        // from a single one)

        labelListList newToOld(invertOneToMany(nNew, oldToNew));

        // Extract points with more than one old one
        pointToMaster.setSize(mesh.nPoints());
        pointToMaster = -1;

        forAll(newToOld, newi)
        {
            const labelList& oldPoints = newToOld[newi];
            if (oldPoints.size() > 1)
            {
                labelList meshPoints
                (
                    labelUIndList(candidates, oldPoints)
                );
                label masteri = min(meshPoints);
                forAll(meshPoints, i)
                {
                    pointToMaster[meshPoints[i]] = masteri;
                }
            }
        }
    }
}


void Foam::snappyLayerDriver::updatePatch
(
    const labelList& patchIDs,
    const mapPolyMesh& map,
    autoPtr<indirectPrimitivePatch>& pp,
    labelList& newToOldPatchPoints
) const
{
    // Update the pp to be consistent with the new mesh

    fvMesh& mesh = meshRefiner_.mesh();

    autoPtr<indirectPrimitivePatch> newPp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );

    // Map from new back to old patch points
    newToOldPatchPoints.setSize(newPp().nPoints());
    newToOldPatchPoints = -1;
    {
        const Map<label>& baseMap = pp().meshPointMap();
        const labelList& newMeshPoints = newPp().meshPoints();

        forAll(newMeshPoints, newPointi)
        {
            const label newMeshPointi = newMeshPoints[newPointi];
            const label oldMeshPointi =
                map.pointMap()[newMeshPointi];
            const auto iter = baseMap.find(oldMeshPointi);
            if (iter.good())
            {
                newToOldPatchPoints[newPointi] = iter();
            }
        }
    }


    // Reset pp. Note: make sure to use std::move - otherwise it
    // will release the pointer before copying and you get
    // memory error. Same if using autoPtr::reset.
    pp = std::move(newPp);

    // Make sure pp has adressing cached
    (void)pp().meshPoints();
    (void)pp().meshPointMap();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyLayerDriver::snappyLayerDriver
(
    meshRefinement& meshRefiner,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const bool dryRun
)
:
    meshRefiner_(meshRefiner),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch),
    dryRun_(dryRun)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::snappyLayerDriver::mergePatchFacesUndo
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const meshRefinement::FaceMergeType mergeType
)
{
    // Clip to 30 degrees. Not helpful!
    //scalar planarAngle = min(30.0, layerParams.featureAngle());
    scalar planarAngle = layerParams.mergePatchFacesAngle();
    scalar minCos = Foam::cos(degToRad(planarAngle));

    scalar concaveCos = Foam::cos(degToRad(layerParams.concaveAngle()));

    Info<< nl
        << "Merging all faces of a cell" << nl
        << "---------------------------" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < " << planarAngle
        << " degrees"
        << " (cos:" << minCos << ')' << nl
        << "    - as long as the resulting face doesn't become concave"
        << " by more than "
        << layerParams.concaveAngle() << " degrees" << nl
        << "      (0=straight, 180=fully concave)" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh));

    labelList duplicateFace(mesh.nFaces(), -1);
    forAll(couples, i)
    {
        const labelPair& cpl = couples[i];
        duplicateFace[cpl[0]] = cpl[1];
        duplicateFace[cpl[1]] = cpl[0];
    }

    label nChanged = meshRefiner_.mergePatchFacesUndo
    (
        minCos,
        concaveCos,
        meshRefiner_.meshedPatches(),
        motionDict,
        duplicateFace,
        mergeType   // How to merge co-planar patch faces
    );

    nChanged += meshRefiner_.mergeEdgesUndo(minCos, motionDict);
}


void Foam::snappyLayerDriver::addLayers
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& patchIDs,
    const label nAllowableErrors,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();

    // faceZones of type internal or baffle (for merging points across)
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner_.getZones(fzTypes);
    }

    // faceZones of type internal (for checking mesh quality across and
    // merging baffles)
    const labelList internalFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::INTERNAL
            )
        )
    );

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    {
        labelList originatingFaceZone;
        meshRefiner_.createZoneBaffles
        (
            identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone
        );

        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing baffled mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }
    }


    // Duplicate points on faceZones of type boundary. Should normally already
    // be done by snapping phase
    {
        autoPtr<mapPolyMesh> map = meshRefiner_.dupNonManifoldBoundaryPoints();
        if (map)
        {
            const labelList& reverseFaceMap = map->reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }
    }



    // Now we have all patches determine the number of layer per patch
    // This will be layerParams.numLayers() except for faceZones where one
    // side does get layers and the other not in which case we want to
    // suppress movement by explicitly setting numLayers 0
    labelList numLayers(layerParams.numLayers());
    {
        labelHashSet layerIDs(patchIDs);
        forAll(mesh.faceZones(), zonei)
        {
            label mpi, spi;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner_.getFaceZoneInfo
            (
                mesh.faceZones()[zonei].name(),
                mpi,
                spi,
                fzType
            );
            if (hasInfo)
            {
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();
                if (layerIDs.found(mpi) && !layerIDs.found(spi))
                {
                    // Only layers on master side. Fix points on slave side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to master patch " << pbm[mpi].name()
                        << " only. Freezing points on slave patch "
                        << pbm[spi].name() << endl;
                    numLayers[spi] = 0;
                }
                else if (!layerIDs.found(mpi) && layerIDs.found(spi))
                {
                    // Only layers on slave side. Fix points on master side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to slave patch " << pbm[spi].name()
                        << " only. Freezing points on master patch "
                        << pbm[mpi].name() << endl;
                    numLayers[mpi] = 0;
                }
            }
        }
    }


    // Duplicate points on faceZones that layers are added to
    labelList pointToMaster;
    dupFaceZonePoints
    (
        meshRefiner_,
        patchIDs,  // patch indices
        numLayers, // number of layers per patch
        baffles,
        pointToMaster
    );


    // Add layers to patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // Now we have
    // - mesh with optional baffles and duplicated points for faceZones
    //   where layers are to be added
    // - pointToMaster    : correspondence for duplicated points
    // - baffles          : list of pairs of faces


    // Calculate 'base' point extrusion
    autoPtr<indirectPrimitivePatch> pp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );
    // Make sure pp has adressing cached before changing mesh later on
    (void)pp().meshPoints();
    (void)pp().meshPointMap();

    // Global face indices engine
    globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges between coupled
    // faces (all edges between same (global)face indices get extruded).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            *pp
        )
    );


    // Point-wise extrusion data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Displacement for all pp.localPoints. Set to large value so does
    // not trigger the minThickness truncation (see syncPatchDisplacement below)
    vectorField basePatchDisp(pp().nPoints(), vector(GREAT, GREAT, GREAT));

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList basePatchNLayers(pp().nPoints(), Zero);

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode> baseExtrudeStatus(pp().nPoints(), EXTRUDE);

    // Ideal number of cells added
    const label nIdealTotAddedCells = setPointNumLayers
    (
        layerParams,

        numLayers,
        patchIDs,
        pp(),
        edgeGlobalFaces,

        basePatchDisp,
        basePatchNLayers,
        baseExtrudeStatus
    );

    // Overall thickness of all layers
    scalarField baseThickness(pp().nPoints());
    // Truncation thickness - when to truncate layers
    scalarIOField baseMinThickness
    (
        IOobject
        (
            "minThickness",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        pp().nPoints()
    );
    // Expansion ratio
    scalarField baseExpansionRatio(pp().nPoints());
    calculateLayerThickness
    (
        pp(),
        patchIDs,
        layerParams,
        meshRefiner_.meshCutter().cellLevel(),
        basePatchNLayers,
        edge0Len,

        baseThickness,
        baseMinThickness,
        baseExpansionRatio
    );


    // Per cell 0 or number of layers in the cell column it is part of
    labelList cellNLayers(mesh.nCells(), Zero);
    // Per face actual overall layer thickness
    scalarField faceRealThickness(mesh.nFaces(), Zero);
    // Per face wanted overall layer thickness
    scalarField faceWantedThickness(mesh.nFaces(), Zero);
    {
        UIndirectList<scalar>(faceWantedThickness, pp().addressing()) =
            avgPointData(pp(), baseThickness);
    }


    // Per patch point the number of layers to add. Is basePatchNLayers
    // for nOuterIter = 1.
    labelList deltaNLayers
    (
        (basePatchNLayers+layerParams.nOuterIter()-1)
       /layerParams.nOuterIter()
    );

    // Per patch point the sum of added layers so far
    labelList nAddedLayers(basePatchNLayers.size(), 0);


    for (label layeri = 0; layeri < layerParams.nOuterIter(); layeri++)
    {
        // Divide layer addition into outer iterations. E.g. if
        // nOutIter = 2, numLayers is 20 for patchA and 1 for patchB
        // this will
        // - iter0:
        //  - add 10 layers to patchA and 1 to patchB
        //  - layers are finalLayerThickness down to the number of layers
        // - iter1 : add 10 layer to patchA and 0 to patchB


        // Exit if nothing to be added
        const label nToAdd = gSum(deltaNLayers);
        Info<< nl
            << "Outer iteration : " << layeri << nl
            << "-------------------" << endl;
        if (debug)
        {
            Info<< "Layers to add in current iteration : " << nToAdd << endl;
        }
        if (nToAdd == 0)
        {
            break;
        }


        // Determine patches for extruded boundary edges. Calculates if any
        // additional processor patches need to be constructed.

        labelList edgePatchID;
        labelList edgeZoneID;
        boolList edgeFlip;
        labelList inflateFaceID;
        determineSidePatches
        (
            meshRefiner_,
            globalFaces,
            edgeGlobalFaces,
            *pp,

            edgePatchID,
            edgeZoneID,
            edgeFlip,
            inflateFaceID
        );


        // Point-wise extrusion data
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Displacement for all pp.localPoints. Set to large value so does
        // not trigger the minThickness truncation (see syncPatchDisplacement
        // below)
        vectorField patchDisp(basePatchDisp);

        // Number of layers for all pp.localPoints. Note: only valid if
        // extrudeStatus = EXTRUDE.
        labelList patchNLayers(deltaNLayers);

        // Whether to add edge for all pp.localPoints.
        List<extrudeMode> extrudeStatus(baseExtrudeStatus);

        // At this point
        //  - patchDisp is either zero or a large value
        //  - patchNLayers is the overall number of layers
        //  - extrudeStatus is EXTRUDE or NOEXTRUDE

        // Determine (wanted) point-wise overall layer thickness and expansion
        // ratio for this deltaNLayers slice of the overall layers
        scalarField sliceThickness(pp().nPoints());
        {
            forAll(baseThickness, pointi)
            {
                sliceThickness[pointi] = layerParameters::layerThickness
                (
                    basePatchNLayers[pointi],   // overall number of layers
                    baseThickness[pointi],      // overall thickness
                    baseExpansionRatio[pointi], // expansion ratio
                    basePatchNLayers[pointi]
                   -nAddedLayers[pointi]
                   -patchNLayers[pointi],       // start index
                    patchNLayers[pointi]        // nLayers to add
                );
            }
        }


        // Current set of topology changes. (changing mesh clears out
        // polyTopoChange)
        polyTopoChange meshMod(mesh.boundaryMesh().size());

        // Shrink mesh, add layers. Returns with any mesh changes in meshMod
        labelList sliceCellNLayers;
        scalarField sliceFaceRealThickness;

        addLayers
        (
            layerParams,
            layerParams.nLayerIter(),

            // Mesh quality
            motionDict,
            layerParams.nRelaxedIter(),
            nAllowableErrors,

            patchIDs,
            internalFaceZones,
            baffles,
            numLayers,
            nIdealTotAddedCells,

            // Patch information
            globalFaces,
            pp(),
            edgeGlobalFaces,
            edgePatchID,
            edgeZoneID,
            edgeFlip,
            inflateFaceID,

            // Per patch point the wanted thickness
            sliceThickness,
            baseMinThickness,
            baseExpansionRatio,

            // Per patch point the wanted&obtained layers and thickness
            patchDisp,
            patchNLayers,
            extrudeStatus,

            // Complete mesh changes
            meshMod,

            // Stats on new mesh
            sliceCellNLayers,       //cellNLayers,
            sliceFaceRealThickness  //faceRealThickness
        );


        // Exit if nothing added
        const label nTotalAdded = gSum(patchNLayers);
        if (debug)
        {
            Info<< nl << "Added in current iteration : " << nTotalAdded
                << " out of : " << gSum(deltaNLayers) << endl;
        }
        if (nTotalAdded == 0)
        {
            break;
        }


        // Update wanted layer statistics
        forAll(patchNLayers, pointi)
        {
            nAddedLayers[pointi] += patchNLayers[pointi];

            if (patchNLayers[pointi] == 0)
            {
                // No layers were added. Make sure that overall extrusion
                // gets reset as well
                unmarkExtrusion
                (
                    pointi,
                    basePatchDisp,
                    basePatchNLayers,
                    baseExtrudeStatus
                );
                basePatchNLayers[pointi] = nAddedLayers[pointi];
                deltaNLayers[pointi] = 0;
            }
            else
            {
                // Adjust the number of layers for the next iteration.
                // Can never be
                // higher than the adjusted overall number of layers.
                // Note: is min necessary?
                deltaNLayers[pointi] = max
                (
                    0,
                    min
                    (
                        deltaNLayers[pointi],
                        basePatchNLayers[pointi] - nAddedLayers[pointi]
                    )
                );
            }
        }


        // At this point we have a (shrunk) mesh and a set of topology changes
        // which will make a valid mesh with layer. Apply these changes to the
        // current mesh.

        {
            // Apply the stored topo changes to the current mesh.
            autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
            mapPolyMesh& map = *mapPtr;

            // Hack to remove meshPhi - mapped incorrectly. TBD.
            mesh.moving(false);
            mesh.clearOut();

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map.hasMotionPoints())
            {
                mesh.movePoints(map.preMotionPoints());
            }
            else
            {
                // Delete mesh volumes.
                mesh.clearOut();
            }

            // Reset the instance for if in overwrite mode
            mesh.setInstance(meshRefiner_.timeName());

            meshRefiner_.updateMesh(map, labelList(0));

            // Update numbering of accumulated cells
            cellNLayers.setSize(map.nOldCells(), 0);
            meshRefinement::updateList
            (
                map.cellMap(),
                label(0),
                cellNLayers
            );
            forAll(cellNLayers, i)
            {
                cellNLayers[i] += sliceCellNLayers[i];
            }

            faceRealThickness.setSize(map.nOldFaces(), scalar(0));
            meshRefinement::updateList
            (
                map.faceMap(),
                scalar(0),
                faceRealThickness
            );
            faceRealThickness += sliceFaceRealThickness;

            meshRefinement::updateList
            (
                map.faceMap(),
                scalar(0),
                faceWantedThickness
            );

            // Print data now that we still have patches for the zones
            //if (meshRefinement::outputLevel() & meshRefinement::OUTPUTLAYERINFO)
            printLayerData
            (
                mesh,
                patchIDs,
                cellNLayers,
                faceWantedThickness,
                faceRealThickness,
                layerParams
            );


            // Dump for debugging
            if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing mesh with layers but disconnected to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    mesh.time().path()/meshRefiner_.timeName()
                );
            }


            // Map baffles, pointToMaster
            mapFaceZonePoints(meshRefiner_, map, baffles, pointToMaster);

            // Map patch and layer settings
            labelList newToOldPatchPoints;
            updatePatch(patchIDs, map, pp, newToOldPatchPoints);

            // Global face indices engine
            globalFaces.reset(mesh.nFaces());

            // Patch-edges to global faces using them
            edgeGlobalFaces = addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                pp()
            );

            // Map patch-point based data
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                vector::uniform(-1),
                basePatchDisp
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                label(0),
                basePatchNLayers
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                extrudeMode::NOEXTRUDE,
                baseExtrudeStatus
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                scalar(0),
                baseThickness
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                scalar(0),
                baseMinThickness
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                GREAT,
                baseExpansionRatio
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                label(0),
                deltaNLayers
            );
            meshRefinement::updateList
            (
                newToOldPatchPoints,
                label(0),
                nAddedLayers
            );
        }
    }


    // Merge baffles
    mergeFaceZonePoints
    (
        pointToMaster,  // per new mesh point : -1 or index of duplicated point
        cellNLayers,        // per new cell : number of layers
        faceRealThickness,  // per new face : actual thickness
        faceWantedThickness // per new face : wanted thickness
    );


    // Do final balancing
    // ~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Balance. No restriction on face zones and baffles.
        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            false,
            false,
            labelList::null(),
            scalarField(mesh.nCells(), 1.0),
            decomposer,
            distributor
        );

        // Re-distribute flag of layer faces and cells
        map().distributeCellData(cellNLayers);
        map().distributeFaceData(faceWantedThickness);
        map().distributeFaceData(faceRealThickness);
    }


    // Write mesh data
    // ~~~~~~~~~~~~~~~

    if (!dryRun_)
    {
        writeLayerData
        (
            mesh,
            patchIDs,
            cellNLayers,
            faceWantedThickness,
            faceRealThickness
        );
    }
}


void Foam::snappyLayerDriver::doLayers
(
    const dictionary& shrinkDict,
    const dictionary& motionDict,
    const layerParameters& layerParams,
    const meshRefinement::FaceMergeType mergeType,
    const bool preBalance,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    addProfiling(layers, "snappyHexMesh::layers");
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Shrinking and layer addition phase" << nl
        << "----------------------------------" << nl
        << endl;


    Info<< "Using mesh parameters " << motionDict << nl << endl;

    // Merge coplanar boundary faces
    if
    (
        mergeType == meshRefinement::FaceMergeType::GEOMETRIC
     || mergeType == meshRefinement::FaceMergeType::IGNOREPATCH
    )
    {
        mergePatchFacesUndo(layerParams, motionDict, mergeType);
    }


    // Per patch the number of layers (-1 or 0 if no layer)
    const labelList& numLayers = layerParams.numLayers();

    // Patches that need to get a layer
    DynamicList<label> patchIDs(numLayers.size());
    label nFacesWithLayers = 0;
    forAll(numLayers, patchi)
    {
        if (numLayers[patchi] > 0)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (!pp.coupled())
            {
                patchIDs.append(patchi);
                nFacesWithLayers += mesh.boundaryMesh()[patchi].size();
            }
            else
            {
                WarningInFunction
                    << "Ignoring layers on coupled patch " << pp.name()
                    << endl;
            }
        }
    }

    // Add contributions from faceZones that get layers
    const faceZoneMesh& fZones = mesh.faceZones();
    forAll(fZones, zonei)
    {
        label mpi, spi;
        surfaceZonesInfo::faceZoneType fzType;
        meshRefiner_.getFaceZoneInfo(fZones[zonei].name(), mpi, spi, fzType);

        if (numLayers[mpi] > 0)
        {
            nFacesWithLayers += fZones[zonei].size();
        }
        if (numLayers[spi] > 0)
        {
            nFacesWithLayers += fZones[zonei].size();
        }
    }


    patchIDs.shrink();

    if (!returnReduceOr(nFacesWithLayers))
    {
        Info<< nl << "No layers to generate ..." << endl;
    }
    else
    {
        // Check that outside of mesh is not multiply connected.
        checkMeshManifold();

        // Check initial mesh
        Info<< "Checking initial mesh ..." << endl;
        labelHashSet wrongFaces(mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces, dryRun_);
        const label nInitErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;


        bool faceZoneOnCoupledFace = false;

        if (!preBalance)
        {
            // Check if there are faceZones on processor boundaries. This
            // requires balancing to move them off the processor boundaries.

            // Is face on a faceZone
            bitSet isExtrudedZoneFace(mesh.nFaces());
            {
                // Add contributions from faceZones that get layers
                const faceZoneMesh& fZones = mesh.faceZones();
                forAll(fZones, zonei)
                {
                    const faceZone& fZone = fZones[zonei];
                    const word& fzName = fZone.name();

                    label mpi, spi;
                    surfaceZonesInfo::faceZoneType fzType;
                    meshRefiner_.getFaceZoneInfo(fzName, mpi, spi, fzType);

                    if (numLayers[mpi] > 0 || numLayers[spi])
                    {
                        isExtrudedZoneFace.set(fZone);
                    }
                }
            }

            bitSet intOrCoupled
            (
                syncTools::getInternalOrCoupledFaces(mesh)
            );

            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                facei++
            )
            {
                if (intOrCoupled[facei] && isExtrudedZoneFace.test(facei))
                {
                    faceZoneOnCoupledFace = true;
                    break;
                }
            }

            Pstream::reduceOr(faceZoneOnCoupledFace);
        }

        // Balance
        if (Pstream::parRun())
        {
            // Calculate wanted number of cells after adding layers
            // (expressed as weight to be passed into decompositionMethod)

            scalarField cellWeights(mesh.nCells(), 1);
            forAll(numLayers, patchi)
            {
                if (numLayers[patchi] > 0)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchi];
                    forAll(pp.faceCells(), i)
                    {
                        cellWeights[pp.faceCells()[i]] += numLayers[patchi];
                    }
                }
            }

            // Add contributions from faceZones that get layers
            const faceZoneMesh& fZones = mesh.faceZones();
            forAll(fZones, zonei)
            {
                const faceZone& fZone = fZones[zonei];
                const word& fzName = fZone.name();

                label mpi, spi;
                surfaceZonesInfo::faceZoneType fzType;
                meshRefiner_.getFaceZoneInfo(fzName, mpi, spi, fzType);

                if (numLayers[mpi] > 0)
                {
                    // Get the owner side for unflipped faces, neighbour side
                    // for flipped ones
                    const labelList& cellIDs = fZone.slaveCells();
                    forAll(cellIDs, i)
                    {
                        if (cellIDs[i] >= 0)
                        {
                            cellWeights[cellIDs[i]] += numLayers[mpi];
                        }
                    }
                }
                if (numLayers[spi] > 0)
                {
                    const labelList& cellIDs = fZone.masterCells();
                    forAll(cellIDs, i)
                    {
                        if (cellIDs[i] >= 0)
                        {
                            cellWeights[cellIDs[i]] += numLayers[mpi];
                        }
                    }
                }
            }


            // Print starting mesh
            Info<< nl;
            meshRefiner_.printMeshInfo
            (
                false,
                "Before layer addition",
                false           //printCellLevel
            );

            // Print expected mesh (if all layers added)
            {
                const scalar nNewCells = sum(cellWeights);
                const scalar nNewCellsAll =
                    returnReduce(nNewCells, sumOp<scalar>());
                const scalar nIdealNewCells = nNewCellsAll / Pstream::nProcs();
                const scalar unbalance = returnReduce
                (
                    mag(1.0-nNewCells/nIdealNewCells),
                    maxOp<scalar>()
                );

                Info<< "Ideal layer addition"
                    << " : cells:" << nNewCellsAll
                    << "  unbalance:" << unbalance
                    << nl << endl;
            }


            if (preBalance || faceZoneOnCoupledFace)
            {
                Info<< nl
                    << "Doing initial balancing" << nl
                    << "-----------------------" << nl
                    << endl;

                // Balance mesh (and meshRefinement). Restrict faceZones to
                // be on internal faces only since they will be converted into
                // baffles.
                autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
                (
                    true,           // keepZoneFaces
                    false,
                    labelList::null(),
                    cellWeights,
                    decomposer,
                    distributor
                );
            }
        }


        // Do all topo changes
        if (layerParams.nOuterIter() == -1)
        {
            // For testing. Can be removed once addLayers below works
            addLayersSinglePass
            (
                layerParams,
                motionDict,
                patchIDs,
                nInitErrors,
                decomposer,
                distributor
            );
        }
        else
        {
            addLayers
            (
                layerParams,
                motionDict,
                patchIDs,
                nInitErrors,
                decomposer,
                distributor
            );
        }
    }
}


// ************************************************************************* //
