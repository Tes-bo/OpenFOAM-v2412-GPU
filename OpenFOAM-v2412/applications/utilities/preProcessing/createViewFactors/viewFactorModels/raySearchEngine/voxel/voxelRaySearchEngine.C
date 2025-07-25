/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "voxelRaySearchEngine.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(voxel, 0);
    addToRunTimeSelectionTable(raySearchEngine, voxel, mesh);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::VF::voxel::setTriangulation(const fvMesh& mesh)
{
    Info<< "\nCreating triangulated surface" << endl;

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(mesh.nBoundaryFaces());
    DynamicList<label> globalFaces(mesh.nBoundaryFaces());
    label nFace = 0;

    const auto& pbm = mesh.boundaryMesh();

    forAll(patchIDs_, i)
    {
        const label patchi = patchIDs_[i];
        const polyPatch& patch = pbm[patchi];
        const pointField& points = patch.points();

        for (const face& f : patch)
        {
            label nTri = 0;
            faceList triFaces(f.nTriangles(points));
            f.triangles(points, nTri, triFaces);

            const label globalFacei =
                globalNumbering_.toGlobal(Pstream::myProcNo(), nFace++);

            for (const face& f : triFaces)
            {
                triangles.push_back(labelledTri(f[0], f[1], f[2], i));
                globalFaces.push_back(globalFacei);
            }
        }
    }

    triToGlobalFace_.transfer(globalFaces);

    Info<< "    Total number of triangles: "
        << returnReduce(triangles.size(), sumOp<label>())
        << endl;

    triangles.shrink();
    surface_ = triSurface(triangles, mesh.points());
    surface_.compactPoints();
}


Foam::labelList Foam::VF::voxel::setFaceVertexHits
(
    const fvMesh& mesh,
    const labelList& patchIDs
)
{
    labelList vertHits(mesh.points().size(), Zero);

    if (mesh.nSolutionD() == 3)
    {
        const auto& pbm = mesh.boundaryMesh();

        // Create a new triangulation based on the surface agglomeration
        for (const label patchI : patchIDs)
        {
            const polyPatch& patch = pbm[patchI];
            for (const face& f : patch)
            {
                for (const label fpi : f)
                {
                    ++vertHits[fpi];
                }
            }
        }

        for (const auto& pp : pbm)
        {
            const labelList& meshPts = pp.meshPoints();

            if (pp.size())
            {
                if (isA<processorPolyPatch>(pp))
                {
                    // Add all processor patch points
                    for (const label pi : meshPts)
                    {
                        ++vertHits[pi];
                    }
                }
                else
                {
                    // Add boundary points

                    const auto& bndyPts = pp.boundaryPoints();

                    for (const label pi : bndyPts)
                    {
                        ++vertHits[meshPts[pi]];
                    }
                }
            }
        }
    }

    return vertHits;
}


void Foam::VF::voxel::setCoarseTriangulation(const fvMesh& mesh)
{
    Info<< "\nCreating triangulated surface" << endl;


    // Filter out fine mesh points along coarse mesh faces


    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(mesh.nBoundaryFaces());
    DynamicList<label> globalFaces(mesh.nBoundaryFaces());
    labelList vertHits = setFaceVertexHits(mesh, patchIDs_);


    // Only simplifying edges for 3D
    const label nVertMin = mesh.nSolutionD() == 3 ? 2 : 0;

    label nInvalid = 0;
    label nFace = 0;

    const auto& pbm = mesh.boundaryMesh();

    for (const label patchi : patchIDs_)
    {
        const polyPatch& patch = pbm[patchi];
        const pointField& points = patch.points();

        for (const face& f : patch)
        {
            DynamicList<label> faceVerts;
            for (const label fpi : f)
            {
                if (vertHits[fpi] > nVertMin)
                {
                    faceVerts.push_back(fpi);
                }
            }

            if (faceVerts.size() < 3)
            {
                ++nInvalid;
                continue;
            }

            label nTri = 0;
            const face reducedFace(faceVerts);
            faceList triFaces(reducedFace.nTriangles(points));
            reducedFace.triangles(points, nTri, triFaces);

            const label globalFacei =
                globalNumbering_.toGlobal(Pstream::myProcNo(), nFace++);

            for (const face& f : triFaces)
            {
                triangles.push_back(labelledTri(f[0], f[1], f[2], patchi));
                globalFaces.push_back(globalFacei);
            }
        }
    }

    triToGlobalFace_.transfer(globalFaces);

    Info<< "    Total number of triangles: "
        << returnReduce(triangles.size(), sumOp<label>())
        << "\n    Number of invalid (removed) triangles: "
        << returnReduce(nInvalid, sumOp<label>())
        << endl;

    triangles.shrink();
    surface_ = triSurface(triangles, mesh.points());
    surface_.compactPoints();
}


void Foam::VF::voxel::broadcast()
{
    // Every processor has whole surface
    const globalIndex globalFaceIdx(globalIndex::gatherOnly{}, surface_.size());
    const globalIndex globalPointIdx
    (
        globalIndex::gatherOnly{},
        surface_.points().size()
    );

    List<labelledTri> globalSurfaceTris(globalFaceIdx.gather(surface_));
    pointField globalSurfacePoints(globalPointIdx.gather(surface_.points()));
    List<label> globalTriToGlobalFace(globalFaceIdx.gather(triToGlobalFace_));


    for (const label proci : globalPointIdx.allProcs())
    {
        const label offset = globalPointIdx.localStart(proci);

        if (offset)
        {
            for
            (
                labelledTri& tri
             :  globalSurfaceTris.slice(globalFaceIdx.range(proci))
            )
            {
                tri[0] += offset;
                tri[1] += offset;
                tri[2] += offset;
            }
        }
    }

    surface_ =
        triSurface
        (
            std::move(globalSurfaceTris),
            std::move(globalSurfacePoints)
        );

    Pstream::broadcast(surface_);

    triToGlobalFace_ = std::move(globalTriToGlobalFace);
    Pstream::broadcast(triToGlobalFace_);
}


Foam::pointHit Foam::VF::voxel::rayTriIntersect
(
    const label trii,
    const point& origin,
    const vector& direction
) const
{
    // Note: origin was made relative to voxel mesh on entry to hit function
    // - need to convert back into global position to be consistent with
    //   triangles for intersection tests

    const auto& ind = surface_[trii];
    const auto& pts = surface_.points();

    // Note: flipped orientation of triangle (into domain) so that we can use
    // visibility check when performing ray-triangle intersections
    const triPointRef tri(pts[ind[0]], pts[ind[2]], pts[ind[1]]);

    static scalar tolerance = 1e-6;

    return
        tri.intersection
        (
            globalPosition(origin),
            direction,
            intersection::VISIBLE,
            tolerance
        );
}


void Foam::VF::voxel::writeBox
(
    OBJstream& os,
    bool lines,
    const boundBox& bb
) const
{
    os.write(treeBoundBox(bb), lines);
}


void Foam::VF::voxel::writeVoxels(const word& fName) const
{
    if (!UPstream::master()) return;

    OBJstream os(fName);
    Info<< "Writing voxels to " << os.name() << endl;

    boundBox bb;
    const bool lines = true;
    for (label i = 0; i < nijk_[0]; ++i)
    {
        for (label j = 0; j < nijk_[1]; ++j)
        {
            for (label k = 0; k < nijk_[2]; ++k)
            {
                bb.min() = voxelMin(i, j, k);
                bb.max() = voxelMax(i, j, k);
                writeBox(os, lines, bb);
            }
        }
    }

    Info<< "- done" << endl;
}


void Foam::VF::voxel::writeTriBoundBoxes(const word& fName) const
{
    if (!UPstream::master()) return;

    OBJstream os(fName);
    Info<< "Writing triangle boundBoxes to " << os.name() << endl;

    const bool lines = true;
    for (const auto& voxeli : objects_)
    {
        for (const label trii : voxeli)
        {
            writeBox(os, lines, objectBbs_[trii]);
        }
    }

    Info<< "- done" << endl;
}


void Foam::VF::voxel::writeHitObjects
(
    const label voxeli,
    const point& origin,
    const vector& dir
) const
{
    OBJstream os("voxel_" + Foam::name(voxeli) + ".obj");

    // Write voxel
    labelVector ijk = this->ijk(voxeli);

    boundBox voxelBb;
    voxelBb.min() = voxelMin(ijk[0], ijk[1], ijk[2]);
    voxelBb.max() = voxelMax(ijk[0], ijk[1], ijk[2]);

    writeBox(os, true, voxelBb);

    for (const auto& trii : objects_[voxeli])
    {
        writeBox(os, true, objectBbs_[trii]);

        const auto& ind = surface_[trii];
        const auto& pts = surface_.points();
        const triPointRef tri(pts[ind[0]], pts[ind[1]], pts[ind[2]]);
        os.write(tri, false);
    }
}


Foam::pointIndexHit Foam::VF::voxel::hitObject
(
    const label voxeli,
    const point& origin,
    const vector& dir,
    const vector& t,
    const scalar minDistance
) const
{
    if (objects_[voxeli].empty()) return pointIndexHit();

    // boundBox rayBb;
    // rayBb.add(origin);
    // rayBb.add(origin + dir*(dir & t));

    label triHiti = -1;
    //rayBb.add(origin + dir);
    //rayBb.inflate(0.01);

    if (debug > 2)
    {
        writeHitObjects(voxeli, origin, dir);
    }

    // Determine all triangles that intersect with ray
    // - only keep nearest hit

    pointHit closestHit;
    for (const auto& trii : objects_[voxeli])
    {
        // Only perform ray/tri intersection if bound boxes intersect
        //if (objectBbs_[trii].overlaps(rayBb))
        {
            pointHit pi = rayTriIntersect(trii, origin, dir);

            if
            (
                pi.hit()
             && (
                    pi.distance() > minDistance
                 && pi.distance() < closestHit.distance()
                )
            )
            {
                triHiti = trii;
                closestHit = pi;
            }
        }
    }

    return pointIndexHit(closestHit, triHiti);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::voxel::voxel(const fvMesh& mesh, const dictionary& dict)
:
    raySearchEngine(mesh, dict),
    bb0_(),
    span0_(Zero),
    nijk_(Zero),
    dxyz_(Zero),
    nRayPerFace_(dict.get<label>("nRayPerFace")),
    nTriPerVoxelMax_(dict.getOrDefault<label>("nTriPerVoxelMax", 50)),
    depthMax_(dict.getOrDefault<label>("depthMax", 5)),
    objects_(),
    objectBbs_()
{
    if (agglomMeshPtr_)
    {
        setCoarseTriangulation(*agglomMeshPtr_);
    }
    else
    {
        setTriangulation(mesh);
    }

    broadcast();

    objectBbs_.resize_nocopy(surface_.size());

    bb0_.add(surface_.points());
    bb0_.inflate(0.01);
    span0_ = bb0_.span();

    {
        scalar maxDx = span0_.x();
        scalar maxDy = span0_.y();
        scalar maxDz = span0_.z();

        const label refDn = 8;
        scalar maxDim = max(maxDx, max(maxDy, maxDz));

        setVoxelDims
        (
            refDn*maxDx/maxDim,
            refDn*maxDy/maxDim,
            refDn*maxDz/maxDim
        );

        objects_.resize_nocopy(nVoxel());
    }

    label depth = 0;
    label trii = 0;
    voxelise(objects_, trii, depth);

    Info<< "\nCreated voxel mesh: " << nijk_ << endl;

    if ((debug > 3) && UPstream::master())
    {
        writeVoxels("voxels.obj");
        writeTriBoundBoxes("triBoundBoxes.obj");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VF::voxel::refineObjects
(
    List<DynamicList<label>>& objects,
    const label triMax
)
{
    refineVoxelDims();

    if (debug > 2) Pout<< "Refining voxels: n=" << nijk_ << endl;

    List<DynamicList<label>> objectsNew(objects.size()*8);

    for (label trii = 0; trii <= triMax; ++trii)
    {
        addBbToVoxels(objectBbs_[trii], trii, objectsNew);
    }

    objects.transfer(objectsNew);
}


Foam::label Foam::VF::voxel::addBbToVoxels
(
    const boundBox& bb,
    const label trii,
    List<DynamicList<label>>& objects
) const
{
    //Pout<< "addBbToVoxels: trii=" << trii << " n=" << nijk_ << endl;

    const point minbb(localPosition(bb.min()));
    const label i0 = max(0, floor(minbb.x()/dxyz_[0]));
    const label j0 = max(0, floor(minbb.y()/dxyz_[1]));
    const label k0 = max(0, floor(minbb.z()/dxyz_[2]));

    const point maxbb(localPosition(bb.max()));
    const label i1 = min(nijk_[0], ceil(maxbb.x()/dxyz_[0]));
    const label j1 = min(nijk_[1], ceil(maxbb.y()/dxyz_[1]));
    const label k1 = min(nijk_[2], ceil(maxbb.z()/dxyz_[2]));

    label nTriMax = 0;

    for (label ii = i0; ii < i1; ++ii)
    {
        for (label jj = j0; jj < j1; ++jj)
        {
            for (label kk = k0; kk < k1; ++kk)
            {
                const label voxeli = this->voxeli(ii, jj, kk);

                objects[voxeli].push_back(trii);
                nTriMax = max(nTriMax, objects[voxeli].size());
            }
        }
    }

    return nTriMax;
}


void Foam::VF::voxel::voxelise
(
    List<DynamicList<label>>& objects,
    const label trii0,
    const label depth
)
{
    if (debug > 2)
    {
        Pout<< "voxelise - start at tri=" << trii0
            << " depth=" << depth
            << endl;
    }

    const auto& points = surface_.points();

    for (label trii = trii0; trii < surface_.size(); ++trii)
    {
        // Triangle bounding box
        boundBox bb(points, surface_[trii]);
        bb.inflate(0.01);
        objectBbs_[trii] = bb;

        const label nVoxelMax = addBbToVoxels(bb, trii, objects);

        // Number of triangles per voxel - if exceed limit, refine voxels...
        if (nVoxelMax > nTriPerVoxelMax_ && depth < depthMax_)
        {
            refineObjects(objects, trii);
            voxelise(objects, trii + 1, depth + 1);
            break;
        }
    }
}


Foam::pointIndexHit Foam::VF::voxel::hit
(
    const label tri0,
    const vector& direction0
) const
{
    if (tri0 > surface_.size()-1)
    {
        FatalErrorInFunction
            << "Index tri0 out of bounds: " << tri0
            << ". Surface size: " << surface_.size()
            << abort(FatalError);
    }

    return hit(surface_.faceCentres()[tri0], direction0);
}


Foam::pointIndexHit Foam::VF::voxel::hit
(
    const point& origin0,
    const vector& direction0
) const
{
    // Initialise return value
    pointIndexHit hitInfo;

    const point origin(localPosition(origin0));

    if (cmptMin(origin) < 0)
    {
        FatalErrorInFunction
            << "Point located outside voxel mesh"
            << " - possible coarsening problem?"
            << abort(FatalError);
    }

    if (magSqr(direction0) < ROOTVSMALL)
    {
        WarningInFunction
            << "Supplied direction has zero size"
            << endl;

        return hitInfo;
    }

    const vector dir(normalised(direction0));

    labelVector ijk(Zero);
    labelVector step(Zero);
    vector tDelta(vector::max);
    vector tMax(vector::max);

    for (direction d = 0; d < 3; ++d)
    {
        ijk[d] = floor(origin[d]/dxyz_[d]);
        step[d] = sign0(dir[d]);
        if (step[d] != 0)
        {
            tDelta[d] = mag(dxyz_[d]/dir[d]);

            scalar voxelMax = (1 + ijk[d] - neg(dir[d]))*dxyz_[d];
            tMax[d] = (voxelMax - origin[d])/dir[d];
        }
    }

    if (debug > 2)
    {
        Info<< "surfBb:" << boundBox(surface_.points())
            << " bb:" << bb0_
            << " origin" << origin0
            << " voxel_origin:" << origin
            << " ijk:" << ijk
            << " step:" << step
            << " dxyz:" << dxyz_
            << " tDelta:" << tDelta
            << " tMax:" << tMax
            << endl;
    }

    auto traverse = [&](const label i)
    {
        ijk[i] += step[i];
        if (outOfBounds(ijk, i)) return false;
        tMax[i] += tDelta[i];
        return true;
    };


    while (true)
    {
        const label voxeli = this->voxeli(ijk);

        if (debug > 2)
        {
            Info<< "ijk:" << ijk
                << " voxeli:" << voxeli
                << " t:" << tMax
                << " objs:" << objects_[voxeli].size()
                << endl;
        }

        hitInfo = hitObject(voxeli, origin, dir, tMax);

        if (hitInfo.hit())
        {
            // Valid hit
            break;
        }
        else
        {
            if (tMax[0] < tMax[1] && tMax[0] < tMax[2])
            {
                if (!traverse(0)) break;
            }
            else if (tMax[1] < tMax[2])
            {
                if (!traverse(1)) break;
            }
            else
            {
                if (!traverse(2)) break;
            }
        }
    }

    return hitInfo;
}


void Foam::VF::voxel::shootRays
(
    labelList& rayStartFaceOut,
    labelList& rayEndFaceOut
) const
{
    (void)mesh_.time().cpuTimeIncrement();

    const pointField& myFc = allCf_[Pstream::myProcNo()];
    const vectorField& myArea = allSf_[Pstream::myProcNo()];

    const label nTotalRays = myFc.size()*nRayPerFace_;
    if (nTotalRays > maxDynListLength)
    {
        FatalErrorInFunction
            << "Dynamic list needs more capacity to support " << nRayPerFace_
            << " rays per face (" << nTotalRays << "). "
            << "Limit set by maxDynListLength: " << maxDynListLength
            << abort(FatalError);
    }

    DynamicList<label> rayStartFace(nTotalRays);
    DynamicList<label> rayEndFace(rayStartFace.size());

    DynamicList<label> endFacei(nTotalRays);
    DynamicList<label> startFacei(nTotalRays);

    const pointField hemi(createHemiPoints(nRayPerFace_));

    Info<< "\nShooting rays" << endl;

    const scalar nDiv = 10;
    const label nStep = floor(myFc.size()/nDiv);
    labelHashSet localFaceHits;

    for (label stepi=0; stepi<nDiv; ++stepi)
    {
        const label step0 = stepi*nStep;
        const label step1 = stepi == nDiv - 1 ? myFc.size() : (stepi+1)*nStep;

        for (label i = step0; i < step1; ++i)
        {
            // Info<< "i/N = " << i << "/" << (myFc.size()-1)
            //     << " step0:" << step0 << " step1:" << step1
            //     << endl;

            localFaceHits.clear();

            const point& origin = myFc[i];
            const vector dir(-normalised(myArea[i]));

            const coordSystem::cartesian cs = createCoordSystem(origin, dir);

            const vectorField pts(cs.transformPoint(hemi));

            for (const auto& p : pts)
            {
                const pointIndexHit hitInfo = hit(origin, p-origin);

                if (hitInfo.hit())
                {
                    label hitFacei = triToGlobalFace_[hitInfo.index()];

                    if (localFaceHits.insert(hitFacei))
                    {
                        endFacei.push_back(hitFacei);
                        startFacei.push_back(i);
                    }
                }
                else
                {
                    // Miss
                }
            }
        }

        const label globalStepi = returnReduce(stepi, minOp<label>()) + 1;

        Info<< "    " << globalStepi/nDiv*100 << "% complete" << endl;
    }


    // Ensure all rays are unique/filter out duplicates
    // - add symmetric connections for non-self-intersecting rays

    if (UPstream::parRun())
    {
        edgeHashSet uniqueRays;
        List<DynamicList<label>> remoteStartFace(Pstream::nProcs());
        List<DynamicList<label>> remoteEndFace(Pstream::nProcs());

        const labelList globalStartFaceIDs
        (
            globalNumbering_.toGlobal(startFacei)
        );

        forAll(globalStartFaceIDs, rayi)
        {
            label start = globalStartFaceIDs[rayi];
            label end = endFacei[rayi];

            const label endProci = globalNumbering_.whichProcID(end);

            if (start > end) Swap(start, end);

            if (uniqueRays.insert(edge(start, end)))
            {
                if (endProci != UPstream::myProcNo())
                {
                    // Convert local face into global face and vice-versa
                    remoteStartFace[endProci].push_back(start);
                    remoteEndFace[endProci].push_back(end);
                }
            }
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send remote data
        for (const int domain : Pstream::allProcs())
        {
            if (domain != Pstream::myProcNo())
            {
                UOPstream str(domain, pBufs);
                str << remoteStartFace[domain]
                    << remoteEndFace[domain];
            }
        }

        pBufs.finishedSends();

        // Consume
        for (const int domain : Pstream::allProcs())
        {
            if (domain != Pstream::myProcNo())
            {
                UIPstream is(domain, pBufs);
                is  >> remoteStartFace[domain]
                    >> remoteEndFace[domain];

                forAll(remoteStartFace[domain], i)
                {
                    const label start = remoteStartFace[domain][i];
                    const label end = remoteEndFace[domain][i];
                    uniqueRays.insert(edge(start, end));
                }
            }
        }

        // Populate ray addressing based on uniqueRays
        for (const edge& e : uniqueRays)
        {
            const label start = e.first();
            const label end = e.second();

            bool sameFace = start == end;

            if (globalNumbering_.isLocal(start))
            {
                // Ray originates from this processor
                const label localStart = globalNumbering_.toLocal(start);

                rayStartFace.append(localStart);
                rayEndFace.append(end);
            }

            if (!sameFace && globalNumbering_.isLocal(end))
            {
                // Ray hits this processor
                // - add symmetric ray from end->start
                const label localEnd = globalNumbering_.toLocal(end);

                rayStartFace.append(localEnd);
                rayEndFace.append(start);
            }
        }
    }
    else
    {
        // Populate ray addressing based on uniqueRays

        edgeHashSet uniqueRays;

        forAll(startFacei, rayi)
        {
            label start = startFacei[rayi];
            label end = endFacei[rayi];

            if (start > end) Swap(start, end);

            if (uniqueRays.insert(edge(start, end)))
            {
                rayStartFace.append(start);
                rayEndFace.append(end);

                if (start != end)
                {
                    // Add symmetric ray from end->start
                    rayStartFace.append(end);
                    rayEndFace.append(start);
                }
            }
        }
    }

    rayStartFaceOut.transfer(rayStartFace);
    rayEndFaceOut.transfer(rayEndFace);

    Info<< "    Time taken: " << mesh_.time().cpuTimeIncrement() << " s"
        << endl;
}


// ************************************************************************* //
