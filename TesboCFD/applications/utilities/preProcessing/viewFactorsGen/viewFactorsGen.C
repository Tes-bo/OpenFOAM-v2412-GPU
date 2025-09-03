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
    viewFactorsGen

Group
    grpPreProcessingUtilities

Description
    This view factors generation application uses a combined approach of
    double area integral (2AI) and double linear integral (2LI). 2AI is used
    when the two surfaces are 'far' apart and 2LI when they are 'close'.
    2LI is integrated along edges using Gaussian quadrature.
    The distance between faces is calculating a ratio between averaged areas
    and the distance between face centres.

    The input from viewFactorsDict are:
    \verbatim
        GaussQuadTol              0.1;      // GaussQuad  error
        distTol                   8;        // R/Average(rm)
        alpha                     0.22;     // Use for common edges for 2LI
    \endverbatim


    For debugging purposes, the following entries can be set in viewFactorsDict:
    \verbatim
        writeViewFactorMatrix     true;
        writeFacesAgglomeration   false;
        dumpRays                  false;

        writeViewFactorMatrix   writes the sum of the VF on each face.
        writeFacesAgglomeration writes the agglomeration
        dumpRays                dumps rays
    \endverbatim

    The participating patches in the VF calculation have to be in the
    'viewFactorWall' patch group (in the polyMesh/boundary file), e.g.

    \verbatim
    floor
    {
        type            wall;
        inGroups        2(wall viewFactorWall);
        nFaces          100;
        startFace       3100;
    }
    \endverbatim

    Compile with -DNO_CGAL only if no CGAL present - CGAL AABB tree performs
    better than the built-in octree.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMesh.H"
#include "meshTools.H"
#include "constants.H"

#include "indirectPrimitivePatch.H"
#include "DynamicField.H"

#include "scalarMatrices.H"
#include "labelListIOList.H"
#include "scalarListIOList.H"

#include "singleCellFvMesh.H"
#include "IOmapDistribute.H"

#ifndef NO_CGAL

// Silence boost bind deprecation warnings (before CGAL-5.2.1)
#include "CGAL/version.h"
#if defined(CGAL_VERSION_NR) && (CGAL_VERSION_NR < 1050211000)
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#endif
#pragma clang diagnostic ignored "-Wbitwise-instead-of-logical"

#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Direction_3 Vector3C;
typedef K::Triangle_3 Triangle;
typedef K::Segment_3 Segment;

typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

// Used boost::optional prior to CGAL-6.0
#if (CGAL_VERSION_NR < 1060000910)
typedef boost::optional
#else
typedef std::optional
#endif
<
    Tree::Intersection_and_primitive_id<Segment>::Type
> Segment_intersection;

#endif // NO_CGAL

using namespace Foam;
using namespace Foam::constant;
using namespace Foam::constant::mathematical;

triSurface triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const labelListIOList& finalAgglom,
    labelList& triSurfaceToAgglom,
    const globalIndex& globalNumbering,
    const polyBoundaryMesh& coarsePatches
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(mesh.nBoundaryFaces());

    label newPatchI = 0;
    label localTriFaceI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            faceList triFaces(f.nTriangles(points));

            label nTri = 0;

            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triFaceI)
            {
                const face& f = triFaces[triFaceI];

                triangles.append(labelledTri(f[0], f[1], f[2], newPatchI));

                nTriTotal++;

                triSurfaceToAgglom[localTriFaceI++] = globalNumbering.toGlobal
                (
                    Pstream::myProcNo(),
                    finalAgglom[patchI][patchFaceI]
                  + coarsePatches[patchI].start()
                );
            }
        }

        newPatchI++;
    }

    //triSurfaceToAgglom.resize(localTriFaceI-1);

    triangles.shrink();
    triSurface surface(triangles, mesh.points());
    surface.compactPoints();


#ifndef NO_CGAL

    // CGAL : every processor has whole surface

    const globalIndex globalFaceIdx
    (
        globalIndex::gatherOnly{},
        surface.size()
    );
    const globalIndex globalPointIdx
    (
        globalIndex::gatherOnly{},
        surface.points().size()
    );

    List<labelledTri> globalSurfaceTris(globalFaceIdx.gather(surface));
    pointField globalSurfacePoints(globalPointIdx.gather(surface.points()));

    //label offset = 0;
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

    surface =
        triSurface
        (
            std::move(globalSurfaceTris),
            std::move(globalSurfacePoints)
        );

    Pstream::broadcast(surface);
#endif

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];

        surface.patches()[newPatchI].index() = patchI;
        surface.patches()[newPatchI].name() = patch.name();
        surface.patches()[newPatchI].geometricType() = patch.type();

        newPatchI++;
    }

    return surface;
}


void writeRays
(
    const fileName& fName,
    const pointField& compactCf,
    const pointField& myFc,
    const labelListList& visibleFaceFaces
)
{
    OFstream str(fName);
    label vertI = 0;

    Pout<< "Dumping rays to " << str.name() << endl;

    forAll(myFc, faceI)
    {
        const labelList visFaces = visibleFaceFaces[faceI];
        forAll(visFaces, faceRemote)
        {
            label compactI = visFaces[faceRemote];
            const point& remoteFc = compactCf[compactI];

            meshTools::writeOBJ(str, myFc[faceI]);
            vertI++;
            meshTools::writeOBJ(str, remoteFc);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }
    str.flush();
}


scalar calculateViewFactorFij2AI
(
    const vector& i,
    const vector& j,
    const vector& dAi,
    const vector& dAj
)
{
    vector r = i - j;
    scalar rMag = mag(r);

    if (rMag > SMALL)
    {
        scalar dAiMag = mag(dAi);
        scalar dAjMag = mag(dAj);

        vector ni = dAi/dAiMag;
        vector nj = dAj/dAjMag;
        scalar cosThetaJ = mag(nj & r)/rMag;
        scalar cosThetaI = mag(ni & r)/rMag;

        return
        (
            (cosThetaI*cosThetaJ*dAjMag*dAiMag)
           /(sqr(rMag)*constant::mathematical::pi)
        );
    }
    else
    {
        return 0;
    }
}


void insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label fromProcI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& matrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(fromProcI, faceI);
        forAll(globalFaces, i)
        {
            matrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


scalar GaussQuad
(

    const scalarList& w,
    const scalarList& p,
    const scalar& magSi,
    const scalar& magSj,
    const vector& di,
    const vector& dj,
    const vector& ci,
    const vector& cj,
    const scalar cosij,
    const scalar alpha,
    label gi
)
{
    scalar dIntFij = 0;
    if (gi == 0)
    {
        vector r(ci - cj);
        if (mag(r) < SMALL)
        {
            r = (alpha*magSi)*di;
        }
        dIntFij = max(cosij*Foam::log(r&r)*magSi*magSj, 0);
    }
    else
    {
        List<vector> pi(w.size());
        forAll (pi, i)
        {
            pi[i] = ci + p[i]*(magSi/2)*di;
        }

        List<vector> pj(w.size());
        forAll (pj, i)
        {
            pj[i] = cj + p[i]*(magSj/2)*dj;
        }

        forAll (w, i)
        {
            forAll (w, j)
            {
                vector r(pi[i] - pj[j]);

                if (mag(r) < SMALL)
                {
                    r = (alpha*magSi)*di;
                     dIntFij +=
                        cosij*w[i]*w[j]*Foam::log(r&r);
                }
                else
                {
                    dIntFij +=
                        cosij*w[i]*w[j]*Foam::log(r&r);
                }

            }
        }

        dIntFij *= (magSi/2) * (magSj/2);

    }
    return dIntFij;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate view factors from face agglomeration array."
        " The finalAgglom generated by faceAgglomerate utility."
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Read view factor dictionary
    IOdictionary viewFactorDict
    (
       IOobject
       (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    const word viewFactorWall("viewFactorWall");

    const bool writeViewFactors =
        viewFactorDict.getOrDefault("writeViewFactorMatrix", false);

    const bool dumpRays =
        viewFactorDict.getOrDefault("dumpRays", false);

    const label debug = viewFactorDict.getOrDefault<label>("debug", 0);

    const scalar GaussQuadTol =
        viewFactorDict.getOrDefault<scalar>("GaussQuadTol", 0.01);

    const scalar distTol =
         viewFactorDict.getOrDefault<scalar>("distTol", 8);

    const scalar alpha =
         viewFactorDict.getOrDefault<scalar>("alpha", 0.21);

    const scalar intTol =
         viewFactorDict.getOrDefault<scalar>("intTol", 1e-2);

    bool useAgglomeration(true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList viewFactorsPatches(patches.indices(viewFactorWall));

    // Read agglomeration map
    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (!finalAgglom.typeHeaderOk<labelListIOList>())
    {
        finalAgglom.setSize(patches.size());
        for (label patchi=0;  patchi < patches.size(); patchi++)
        {
            finalAgglom[patchi] = identity(patches[patchi].size());
        }
        useAgglomeration = false;
    }

    // Create the coarse mesh  using agglomeration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout << "\nCreating single cell mesh..." << endl;
    }

    singleCellFvMesh coarseMesh
    (
        IOobject
        (
            "coarse:" + mesh.name(),
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        finalAgglom
    );

    if (debug)
    {
        Pout << "\nCreated single cell mesh..." << endl;
    }


    // Calculate total number of fine and coarse faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nCoarseFaces = 0;      //total number of coarse faces
    label nFineFaces = 0;        //total number of fine faces

    const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

    for (const label patchi : viewFactorsPatches)
    {
        nCoarseFaces += coarsePatches[patchi].size();
        nFineFaces += patches[patchi].size();
    }


    Info<< "\nTotal number of coarse faces: "
        << returnReduce(nCoarseFaces, sumOp<label>())
        << endl;

    if (Pstream::master() && debug)
    {
        Pout << "\nView factor patches included in the calculation : "
             << viewFactorsPatches << endl;
    }

    // Collect local Cf and Sf on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<point> localCoarseCf(nCoarseFaces);
    DynamicList<point> localCoarseSf(nCoarseFaces);
    DynamicList<label> localAgg(nCoarseFaces);
    labelHashSet includePatches;

    for (const label patchID : viewFactorsPatches)
    {
        const polyPatch& pp = patches[patchID];
        const labelList& agglom = finalAgglom[patchID];

        includePatches.insert(patchID);

        if (agglom.size() > 0)
        {
            label nAgglom = max(agglom)+1;
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchID];

            const pointField& coarseCf =
                coarseMesh.Cf().boundaryField()[patchID];
            const pointField& coarseSf =
                coarseMesh.Sf().boundaryField()[patchID];

            forAll(coarseCf, faceI)
            {
                point cf = coarseCf[faceI];

                const label coarseFaceI = coarsePatchFace[faceI];
                const labelList& fineFaces = coarseToFine[coarseFaceI];
                const label agglomI =
                    agglom[fineFaces[0]] + coarsePatches[patchID].start();

                // Construct single face
                uindirectPrimitivePatch upp
                (
                    UIndirectList<face>(pp, fineFaces),
                    pp.points()
                );

                List<point> availablePoints
                (
                    upp.faceCentres().size()
                  + upp.localPoints().size()
                );

                SubList<point>
                (
                    availablePoints,
                    upp.faceCentres().size()
                ) = upp.faceCentres();

                SubList<point>
                (
                    availablePoints,
                    upp.localPoints().size(),
                    upp.faceCentres().size()
                ) = upp.localPoints();

                point cfo = cf;
                scalar dist = GREAT;
                forAll(availablePoints, iPoint)
                {
                    point cfFine = availablePoints[iPoint];
                    if (mag(cfFine-cfo) < dist)
                    {
                        dist = mag(cfFine-cfo);
                        cf = cfFine;
                    }
                }

                point sf = coarseSf[faceI];
                localCoarseCf.append(cf);
                localCoarseSf.append(sf);
                localAgg.append(agglomI);

            }
        }
    }

    // Distribute local coarse Cf and Sf for shooting rays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointField> remoteCoarseCf(Pstream::nProcs());
    List<pointField> remoteCoarseSf(Pstream::nProcs());
    List<labelField> remoteCoarseAgg(Pstream::nProcs());

    remoteCoarseCf[Pstream::myProcNo()] = localCoarseCf;
    remoteCoarseSf[Pstream::myProcNo()] = localCoarseSf;
    remoteCoarseAgg[Pstream::myProcNo()] = localAgg;

    Pstream::allGatherList(remoteCoarseCf);
    Pstream::allGatherList(remoteCoarseSf);
    Pstream::allGatherList(remoteCoarseAgg);


    globalIndex globalNumbering(nCoarseFaces);

    // Set up searching engine for obstacles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef NO_CGAL
    // Using octree
    #include "searchingEngine.H"
#else
    // Using CGAL aabbtree (faster, more robust)
    #include "searchingEngine_CGAL.H"
#endif

    // Determine rays between coarse face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    DynamicList<label> rayEndFace(rayStartFace.size());

    // Return rayStartFace in local index and rayEndFace in global index
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef NO_CGAL
    // Using octree, distributedTriSurfaceMesh
    #include "shootRays.H"
#else
    // Using CGAL aabbtree (faster, more robust)
    #include "shootRays_CGAL.H"
#endif

    // Calculate number of visible faces from local index
    labelList nVisibleFaceFaces(nCoarseFaces, Zero);

    forAll(rayStartFace, i)
    {
        nVisibleFaceFaces[rayStartFace[i]]++;
    }

    labelListList visibleFaceFaces(nCoarseFaces);

    label nViewFactors = 0;
    forAll(nVisibleFaceFaces, faceI)
    {
        visibleFaceFaces[faceI].setSize(nVisibleFaceFaces[faceI]);
        nViewFactors += nVisibleFaceFaces[faceI];
    }

    // - Construct compact numbering
    // - return map from remote to compact indices
    //   (per processor (!= myProcNo) a map from remote index to compact index)
    // - construct distribute map
    // - renumber rayEndFace into compact addressing

    List<Map<label>> compactMap(Pstream::nProcs());

    mapDistribute map(globalNumbering, rayEndFace, compactMap);

    // visibleFaceFaces has:
    //    (local face, local viewed face) = compact viewed face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nVisibleFaceFaces = 0;
    forAll(rayStartFace, i)
    {
        label faceI = rayStartFace[i];
        label compactI = rayEndFace[i];
        visibleFaceFaces[faceI][nVisibleFaceFaces[faceI]++] = compactI;
    }

    // Construct data in compact addressing
    // (2AA) need coarse (Ai), fine Sf (dAi) and fine Cf(r) to calculate Fij
    // (2LI) need edges (li)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField compactCoarseCf(map.constructSize(), Zero);
    pointField compactCoarseSf(map.constructSize(), Zero);
    List<List<point>> compactFineSf(map.constructSize());
    List<List<point>> compactFineCf(map.constructSize());

    DynamicList<List<point>> compactPoints(map.constructSize());

    DynamicList<label> compactPatchId(map.constructSize());

    // Insert my coarse local values
    SubList<point>(compactCoarseSf, nCoarseFaces) = localCoarseSf;
    SubList<point>(compactCoarseCf, nCoarseFaces) = localCoarseCf;

    const faceList& faces = mesh.faces();

    // Insert my fine local values
    label compactI = 0;
    forAll(viewFactorsPatches, i)
    {
        label patchID = viewFactorsPatches[i];

        const labelList& agglom = finalAgglom[patchID];
        if (agglom.size() > 0)
        {
            label nAgglom = max(agglom)+1;
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchID];

            const polyPatch& pp = patches[patchID];

            forAll(coarseToFine, coarseI)
            {
                compactPatchId.append(patchID);
                List<point>& fineCf = compactFineCf[compactI];
                List<point>& fineSf = compactFineSf[compactI];

                label startFace = pp.start();

                const vectorField locPoints
                (
                    mesh.points(),
                    faces[coarseI + startFace]
                );

                const label coarseFaceI = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceI];

                fineCf.setSize(fineFaces.size());
                fineSf.setSize(fineFaces.size());

                compactPoints.append(locPoints);

                fineCf = UIndirectList<point>
                (
                    mesh.Cf().boundaryField()[patchID],
                    coarseToFine[coarseFaceI]
                );
                fineSf = UIndirectList<point>
                (
                    mesh.Sf().boundaryField()[patchID],
                    coarseToFine[coarseFaceI]
                );

                compactI++;
            }
        }
    }

    if (Pstream::master() && debug)
    {
        Info<< "map distribute..."  << endl;
    }

    // Do all swapping
    map.distribute(compactCoarseSf);
    map.distribute(compactCoarseCf);
    map.distribute(compactFineCf);
    map.distribute(compactFineSf);
    map.distribute(compactPoints);
    map.distribute(compactPatchId);

    // Plot all rays between visible faces.
    if (dumpRays)
    {
        writeRays
        (
            runTime.path()/"allVisibleFaces.obj",
            compactCoarseCf,
            remoteCoarseCf[Pstream::myProcNo()],
            visibleFaceFaces
        );
    }


    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarListIOList F2LI
    (
        IOobject
        (
            "F",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        nCoarseFaces
    );

    const label totalPatches =
        returnReduce(coarsePatches.size(), maxOp<label>());

    // Matrix sum in j(Fij) for each i (if enclosure sum = 1)
    scalarSquareMatrix sumViewFactorPatch(totalPatches, Zero);

    scalarList patchArea(totalPatches, Zero);

    if (Pstream::master())
    {
        Info<< "\nCalculating view factors..." << endl;
    }

    FixedList<scalarList, 5> GaussPoints;
    GaussPoints[0].setSize(1);
    GaussPoints[0] = 0;

    GaussPoints[1].setSize(2);
    GaussPoints[1][0] =  1/std::sqrt(3);
    GaussPoints[1][1] = -1/std::sqrt(3);

    GaussPoints[2].setSize(3);
    GaussPoints[2][0] =  0;
    GaussPoints[2][1] =  std::sqrt(3.0/5.0);
    GaussPoints[2][2] = -std::sqrt(3.0/5.0);

    GaussPoints[3].setSize(4);
    GaussPoints[3][0] = std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0));
    GaussPoints[3][1] = -GaussPoints[3][0];
    GaussPoints[3][2] = std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0));
    GaussPoints[3][3] = -GaussPoints[3][2];

    GaussPoints[4].setSize(5);
    GaussPoints[4][0] =  0;
    GaussPoints[4][1] = (1.0/3.0)*std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0));
    GaussPoints[4][2] = -GaussPoints[4][1];
    GaussPoints[4][3] = (1.0/3.0)*std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0));
    GaussPoints[4][4] = -GaussPoints[4][3];


    FixedList<scalarList, 5> GaussWeights;
    GaussWeights[0].setSize(1);
    GaussWeights[0] = 2;

    GaussWeights[1].setSize(2);
    GaussWeights[1][0] =  1;
    GaussWeights[1][1] =  1;

    GaussWeights[2].setSize(3);
    GaussWeights[2][0] =  8.0/9.0;
    GaussWeights[2][1] =  5.0/9.0;
    GaussWeights[2][2] =  5.0/9.0;

    GaussWeights[3].setSize(4);
    GaussWeights[3][0] =  (18.0 + std::sqrt(30))/36.0;
    GaussWeights[3][1] =  (18.0 + std::sqrt(30))/36.0;
    GaussWeights[3][2] =  (18.0 - std::sqrt(30))/36.0;
    GaussWeights[3][3] =  (18.0 - std::sqrt(30))/36.0;

    GaussWeights[4].setSize(5);
    GaussWeights[4][0] =  128.0/225.0;
    GaussWeights[4][1] =  (322.0 + 13.0*std::sqrt(70))/900.0;
    GaussWeights[4][2] =  (322.0 + 13.0*std::sqrt(70))/900.0;
    GaussWeights[4][3] =  (322.0 - 13.0*std::sqrt(70))/900.0;
    GaussWeights[4][4] =  (322.0 - 13.0*std::sqrt(70))/900.0;

    const label maxQuadOrder = 5;

    if (mesh.nSolutionD() == 3)
    {
        forAll(localCoarseSf, coarseFaceI)
        {
            const List<point>& localFineSf = compactFineSf[coarseFaceI];
            const vector Ai = sum(localFineSf);
            const List<point>& localFineCf = compactFineCf[coarseFaceI];
            const label fromPatchId = compactPatchId[coarseFaceI];

            const List<point>& lPoints = compactPoints[coarseFaceI];

            patchArea[fromPatchId] += mag(Ai);

            const labelList& visCoarseFaces = visibleFaceFaces[coarseFaceI];

            forAll(visCoarseFaces, visCoarseFaceI)
            {
                //F2AI[coarseFaceI].setSize(visCoarseFaces.size());
                F2LI[coarseFaceI].setSize(visCoarseFaces.size());
                label compactJ = visCoarseFaces[visCoarseFaceI];
                const List<point>& remoteFineSj = compactFineSf[compactJ];
                const List<point>& remoteFineCj = compactFineCf[compactJ];

                const List<point>& rPointsCj = compactPoints[compactJ];

                const label toPatchId = compactPatchId[compactJ];

                bool far(false);
                // Relative distance
                forAll(localFineSf, i)
                {
                    const scalar dAi =
                        Foam::sqrt
                        (
                            mag(localFineSf[i])/constant::mathematical::pi
                        );
                    const vector& dCi = localFineCf[i];

                    forAll(remoteFineSj, j)
                    {
                        const scalar dAj =
                            Foam::sqrt
                            (
                                mag(remoteFineSj[j])/constant::mathematical::pi
                            );
                        const vector& dCj = remoteFineCj[j];

                        const scalar dist = mag(dCi - dCj)/((dAi + dAj)/2);

                        if (dist > distTol)
                        {
                            far = true;
                        }
                    }
                }

                if (far)
                {
                    // 2AI method
                    scalar F2AIij = 0;

                    forAll(localFineSf, i)
                    {
                        const vector& dAi = localFineSf[i];
                        const vector& dCi = localFineCf[i];

                        forAll(remoteFineSj, j)
                        {
                            const vector& dAj = remoteFineSj[j];
                            const vector& dCj = remoteFineCj[j];

                            scalar dIntFij = calculateViewFactorFij2AI
                            (
                                dCi,
                                dCj,
                                dAi,
                                dAj
                            );

                            F2AIij += dIntFij;
                        }
                    }
                    F2LI[coarseFaceI][visCoarseFaceI] = F2AIij/mag(Ai);
                }
                else
                {
                    // 2LI method
                    label nLocal = lPoints.size();
                    label nRemote = rPointsCj.size();

                    // Using sub-divisions (quadrature)
                    scalar oldEToeInt = 0;
                    for (label gi=0; gi < maxQuadOrder; gi++)
                    {
                        scalar F2LIij = 0;
                        for(label i=0; i<nLocal; i++)
                        {
                            vector si;
                            vector ci;

                            vector sj;
                            vector cj;

                            if (i == 0)
                            {
                                si = lPoints[i] - lPoints[nLocal-1];
                                ci = (lPoints[i] + lPoints[nLocal-1])/2;
                            }
                            else
                            {
                                si = lPoints[i] - lPoints[i-1];
                                ci = (lPoints[i] + lPoints[i-1])/2;
                            }

                            for(label j=0; j<nRemote; j++)
                            {
                                if (j == 0)
                                {
                                    sj = rPointsCj[j]-rPointsCj[nRemote-1];
                                    cj = (rPointsCj[j]+rPointsCj[nRemote-1])/2;
                                }
                                else
                                {
                                    sj = rPointsCj[j] - rPointsCj[j-1];
                                    cj = (rPointsCj[j] + rPointsCj[j-1])/2;
                                }


                                scalar magSi = mag(si);
                                scalar magSj = mag(sj);
                                scalar cosij = (si & sj)/(magSi * magSj);

                                vector di = si/magSi;
                                vector dj = sj/magSj;

                                label quadOrder = gi;
                                const vector r(ci - cj);
                                // Common edges use n = 0
                                if (mag(r) < SMALL)
                                {
                                    quadOrder = 0;
                                }

                                scalar dIntFij =
                                    GaussQuad
                                    (
                                        GaussWeights[gi],
                                        GaussPoints[gi],
                                        magSi,
                                        magSj,
                                        di,
                                        dj,
                                        ci,
                                        cj,
                                        cosij,
                                        alpha,
                                        quadOrder
                                    );

                                F2LIij += dIntFij;
                            }
                        }

                        const scalar err
                        (
                            mag(F2LIij) > ROOTVSMALL
                          ? (F2LIij-oldEToeInt)/F2LIij
                          : 0
                        );

                        if
                        (
                            (mag(err) < GaussQuadTol && gi > 0)
                         || gi == maxQuadOrder-1
                        )
                        {
                            F2LI[coarseFaceI][visCoarseFaceI] =
                                F2LIij/mag(Ai)/4/constant::mathematical::pi;
                            break;
                        }
                        else
                        {
                            oldEToeInt = F2LIij;
                        }
                    }
                }

                sumViewFactorPatch[fromPatchId][toPatchId] +=
                    F2LI[coarseFaceI][visCoarseFaceI]*mag(Ai);
            }
        }
    }
    else if (mesh.nSolutionD() == 2)
    {
        const boundBox& box = mesh.bounds();
        const Vector<label>& dirs = mesh.geometricD();
        vector emptyDir = Zero;
        forAll(dirs, i)
        {
            if (dirs[i] == -1)
            {
                emptyDir[i] = 1.0;
            }
        }

        scalar wideBy2 = (box.span() & emptyDir)*2.0;

        forAll(localCoarseSf, coarseFaceI)
        {
            const vector& Ai = localCoarseSf[coarseFaceI];
            const vector& Ci = localCoarseCf[coarseFaceI];
            vector Ain = Ai/mag(Ai);
            vector R1i = Ci + (mag(Ai)/wideBy2)*(Ain ^ emptyDir);
            vector R2i = Ci - (mag(Ai)/wideBy2)*(Ain ^ emptyDir) ;

            const label fromPatchId = compactPatchId[coarseFaceI];
            patchArea[fromPatchId] += mag(Ai);

            const labelList& visCoarseFaces = visibleFaceFaces[coarseFaceI];
            forAll(visCoarseFaces, visCoarseFaceI)
            {
                F2LI[coarseFaceI].setSize(visCoarseFaces.size());
                label compactJ = visCoarseFaces[visCoarseFaceI];
                const vector& Aj = compactCoarseSf[compactJ];
                const vector& Cj = compactCoarseCf[compactJ];

                const label toPatchId = compactPatchId[compactJ];

                vector Ajn = Aj/mag(Aj);
                vector R1j = Cj + (mag(Aj)/wideBy2)*(Ajn ^ emptyDir);
                vector R2j = Cj - (mag(Aj)/wideBy2)*(Ajn ^ emptyDir);

                scalar d1 = mag(R1i - R2j);
                scalar d2 = mag(R2i - R1j);
                scalar s1 = mag(R1i - R1j);
                scalar s2 = mag(R2i - R2j);

                scalar Fij = mag((d1 + d2) - (s1 + s2))/(4.0*mag(Ai)/wideBy2);

                F2LI[coarseFaceI][visCoarseFaceI] = Fij;
                sumViewFactorPatch[fromPatchId][toPatchId] += Fij*mag(Ai);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << " View factors are not available in 1D "
            << exit(FatalError);
    }

    // Write view factors matrix in listlist form
    F2LI.write();

    reduce(sumViewFactorPatch, sumOp<scalarSquareMatrix>());
    reduce(patchArea, sumOp<scalarList>());

    if (Pstream::master() && debug)
    {
        forAll(viewFactorsPatches, i)
        {
            label patchI =  viewFactorsPatches[i];
            for (label j=i; j<viewFactorsPatches.size(); j++)
            {
                label patchJ =  viewFactorsPatches[j];

                Info << "F" << patchI << patchJ << ": "
                     << sumViewFactorPatch[patchI][patchJ]/patchArea[patchI]
                     << endl;
            }
        }
    }

    if (writeViewFactors)
    {
        if (Pstream::master())
        {
            Info << "Writing view factor matrix..." << endl;
        }

        volScalarField viewFactorField
        (
            IOobject
            (
                "viewFactorField",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        label compactI = 0;

        volScalarField::Boundary& vfbf = viewFactorField.boundaryFieldRef();
        forAll(viewFactorsPatches, i)
        {
            label patchID = viewFactorsPatches[i];
            const labelList& agglom = finalAgglom[patchID];
            if (agglom.size() > 0)
            {
                label nAgglom = max(agglom)+1;
                labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
                const labelList& coarsePatchFace =
                    coarseMesh.patchFaceMap()[patchID];

                forAll(coarseToFine, coarseI)
                {
                    const scalar FiSum = sum(F2LI[compactI]);

                    const label coarseFaceID = coarsePatchFace[coarseI];
                    const labelList& fineFaces = coarseToFine[coarseFaceID];
                    forAll(fineFaces, fineId)
                    {
                        const label faceID = fineFaces[fineId];
                        vfbf[patchID][faceID] = FiSum;
                    }
                    compactI++;
                }
            }
        }
        viewFactorField.write();
    }


    // Invert compactMap (from processor+localface to compact) to go
    // from compact to processor+localface (expressed as a globalIndex)
    // globalIndex globalCoarFaceNum(coarseMesh.nFaces());
    labelList compactToGlobal(map.constructSize());

    // Local indices first (note: are not in compactMap)
    for (label i = 0; i < globalNumbering.localSize(); i++)
    {
        compactToGlobal[i] = globalNumbering.toGlobal(i);
    }


    forAll(compactMap, procI)
    {
        const Map<label>& localToCompactMap = compactMap[procI];

        forAllConstIters(localToCompactMap, iter)
        {
            compactToGlobal[*iter] = globalNumbering.toGlobal
            (
                procI,
                iter.key()
            );
        }
    }


    labelListList globalFaceFaces(visibleFaceFaces.size());

    // Create globalFaceFaces needed to insert view factors
    // in F to the global matrix Fmatrix
    forAll(globalFaceFaces, faceI)
    {
        globalFaceFaces[faceI] = renumber
        (
            compactToGlobal,
            visibleFaceFaces[faceI]
        );
    }

    labelListIOList IOglobalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        std::move(globalFaceFaces)
    );
    IOglobalFaceFaces.write();


    IOmapDistribute IOmapDist
    (
        IOobject
        (
            "mapDist",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(map)
    );

    IOmapDist.write();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
