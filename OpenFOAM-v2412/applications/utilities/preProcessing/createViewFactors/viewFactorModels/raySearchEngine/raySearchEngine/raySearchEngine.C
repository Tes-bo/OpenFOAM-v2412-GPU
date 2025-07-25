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

#include "raySearchEngine.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "meshTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(raySearchEngine, 0);
    defineRunTimeSelectionTable(raySearchEngine, mesh);

    const label raySearchEngine::maxDynListLength = 1000000000;
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::VF::raySearchEngine::check
(
    const labelList& nVisibleFaceFaces
)
{
    label nRay = 0;
    label nFaceMin = labelMax;
    label nFaceMax = labelMin;
    for (const label n : nVisibleFaceFaces)
    {
        nFaceMin = min(nFaceMin, n);
        nFaceMax = max(nFaceMax, n);
        nRay += n;
    }

    const label nFace = nVisibleFaceFaces.size();
    const label nGlobalRays = returnReduce(nRay, sumOp<label>());

    if (nGlobalRays == 0)
    {
        FatalErrorInFunction
            << "No rays identified - view factors will not be calculated"
            << exit(FatalError);
    }

    const label globalNFacesMin = returnReduce(nFaceMin, minOp<label>());
    const label globalNFacesMax = returnReduce(nFaceMax, maxOp<label>());
    const label nGlobalFaces = returnReduce(nFace, sumOp<label>());
    const scalar avgFaces = nGlobalRays/scalar(nGlobalFaces);

    Info<< "\nRay summary:" << nl
        << "    Number of rays: " << nGlobalRays << nl
        << "    Number of rays-per-face (min, max, average): ("
        << globalNFacesMin << ", "
        << globalNFacesMax << ", "
        << avgFaces << ")" << endl;
}


Foam::label Foam::VF::raySearchEngine::closestPointIndex
(
    const point& p0,
    const List<point>& pts
)
{
    label pointi = -1;

    scalar distSqr = GREAT;
    forAll(pts, pti)
    {
        const scalar d2 = magSqr(pts[pti] - p0);
        if (d2 < distSqr)
        {
            pointi = pti;
            distSqr = d2;
        }
    }

    return pointi;
}


void Foam::VF::raySearchEngine::createAgglomeration(const IOobject& io)
{
    Info<< "\nFace agglomeration: active" << nl
        << "    Reading file " << io.name() << endl;

    // Read agglomeration map
    const labelListIOList finalAgglom(io);

    Info<< "    Creating coarse mesh" << nl;

    agglomMeshPtr_.reset
    (
        new singleCellFvMesh
        (
            IOobject
            (
                IOobject::scopedName("agglom", mesh_.name()),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            finalAgglom
        )
    );

    const auto& coarseMesh = agglomMeshPtr_();


    // Calculate total number of fine and coarse faces

    nCoarseFace_ = 0;
    nFace_ = 0;

    const polyBoundaryMesh& finePatches = mesh_.boundaryMesh();
    const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

    for (const label patchi : patchIDs_)
    {
        nCoarseFace_ += coarsePatches[patchi].size();
        nFace_ += finePatches[patchi].size();
    }

    Info<< "\nTotal number of coarse faces: "
        << returnReduce(nCoarseFace_, sumOp<label>())
        << endl;

    Info<< "\nTotal number of fine faces: "
        << returnReduce(nFace_, sumOp<label>())
        << endl;

    // Collect local Cf, Sf, agglom index on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<point> localCf(nCoarseFace_);
    DynamicList<vector> localSf(nCoarseFace_);
    DynamicList<label> localAgg(nCoarseFace_);

    for (const label patchi : patchIDs_)
    {
        const labelList& agglom = finalAgglom[patchi];

        if (agglom.empty()) continue;

        label nAgglom = max(agglom) + 1;
        const labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
        const labelList& coarsePatchFace = coarseMesh.patchFaceMap()[patchi];

        const pointField& coarseCf = coarseMesh.Cf().boundaryField()[patchi];
        const vectorField& coarseSf = coarseMesh.Sf().boundaryField()[patchi];

        const polyPatch& pp = finePatches[patchi];
        patchAreas_[patchi] += sum(coarseMesh.magSf().boundaryField()[patchi]);

        forAll(coarseCf, facei)
        {
            const label coarseFacei = coarsePatchFace[facei];
            const label agglomi = coarseFacei + coarsePatches[patchi].start();

            // Construct single coarse face
            const labelList& fineFaces = coarseToFine[coarseFacei];
            uindirectPrimitivePatch cf
            (
                UIndirectList<face>(pp, fineFaces),
                pp.points()
            );

            // Collect all points (vertices, face centres)
            const label nFaces = cf.faceCentres().size();
            const label nPoints = cf.localPoints().size();
            List<point> allPoints(nFaces + nPoints);
            SubList<point>(allPoints, nFaces) = cf.faceCentres();
            SubList<point>(allPoints, nPoints, nFaces) = cf.localPoints();

            // Coarse face centre set to closest point
            const label pti = closestPointIndex(coarseCf[facei], allPoints);

            if (pti != -1)
            {
                localCf.push_back(allPoints[pti]);
                localSf.push_back(coarseSf[facei]);
                localAgg.push_back(agglomi);
            }
        }
    }

    Info<< "\nAssembled coarse patch data" << endl;

    // Distribute local coarse Cf and Sf for shooting rays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    allCf_[Pstream::myProcNo()].transfer(localCf);
    allSf_[Pstream::myProcNo()].transfer(localSf);
    allAgg_[Pstream::myProcNo()].transfer(localAgg);

    Pstream::allGatherList(allCf_);
    Pstream::allGatherList(allSf_);
    Pstream::allGatherList(allAgg_);

    Pstream::listCombineGather(patchAreas_, plusEqOp<scalar>());
    Pstream::broadcast(patchAreas_);

    globalNumbering_ = globalIndex(nCoarseFace_);
}


void Foam::VF::raySearchEngine::createGeometry()
{
    DynamicList<point> localCf(mesh_.nBoundaryFaces());
    DynamicList<vector> localSf(mesh_.nBoundaryFaces());

    const auto& pbm = mesh_.boundaryMesh();

    for (const label patchi : patchIDs_)
    {
        localCf.push_back(pbm[patchi].faceCentres());
        localSf.push_back(pbm[patchi].faceAreas());

        patchAreas_[patchi] += sum(mesh_.magSf().boundaryField()[patchi]);
    }

    Info<< "\nAssembled patch data" << endl;

    nFace_ = localCf.size();
    nCoarseFace_ = -1;

    allCf_[Pstream::myProcNo()].transfer(localCf);
    allSf_[Pstream::myProcNo()].transfer(localSf);

    Pstream::allGatherList(allCf_);
    Pstream::allGatherList(allSf_);

    Pstream::listCombineGather(patchAreas_, plusEqOp<scalar>());
    Pstream::broadcast(patchAreas_);

    globalNumbering_ = globalIndex(nFace_);
}


void Foam::VF::raySearchEngine::createParallelAddressing
(
    labelList& rayEndFace
) const
{
    // Construct compact numbering
    // - return map from remote to compact indices
    //   (per processor (!= myProcNo) a map from remote index to compact index)
    // - construct distribute map
    // - renumber rayEndFace into compact addressing

    DebugInfo << "\nCreating map distribute" << endl;

    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset(new mapDistribute(globalNumbering_, rayEndFace, compactMap));

    DebugInfo << "\nCreating compact-to-global addressing" << endl;

    // Invert compactMap (from processor+localface to compact) to go
    // from compact to processor+localface (expressed as a globalIndex)
    compactToGlobal_.resize_nocopy(mapPtr_->constructSize());

    // Local indices first
    // Note: are not in compactMap
    for (label i = 0; i < globalNumbering_.localSize(); ++i)
    {
        compactToGlobal_[i] = globalNumbering_.toGlobal(i);
    }

    forAll(compactMap, proci)
    {
        const Map<label>& localToCompactMap = compactMap[proci];

        forAllConstIters(localToCompactMap, iter)
        {
            compactToGlobal_[*iter] =
                globalNumbering_.toGlobal(proci, iter.key());
        }
    }
}


Foam::coordSystem::cartesian Foam::VF::raySearchEngine::createCoordSystem
(
    const point& origin,
    const vector& dir
) const
{
    vector axis(Zero);

    for (direction d=0; d<3; ++d)
    {
        axis = dir^tensor::I.col(d);

        // Remove empty direction for 2D
        if (mesh_.nSolutionD() == 2)
        {
            meshTools::constrainDirection(mesh_, mesh_.solutionD(), axis);
        }

        if (magSqr(axis) > 0)
        {
            axis.normalise();
            break;
        }
    }

    return coordSystem::cartesian(origin, dir, axis);
}


Foam::tmp<Foam::pointField> Foam::VF::raySearchEngine::createHemiPoints
(
    const label nRayPerFace
) const
{
    auto themiPts = tmp<pointField>::New(nRayPerFace);
    auto& hemiPts = themiPts.ref();

    const label nPoints = hemiPts.size();

    if (mesh_.nSolutionD() == 3)
    {
        // Point in range -1 < x < 1; -1 < y < 1; 0 < z 1

        forAll(hemiPts, pointi)
        {
            const scalar phi = Foam::acos(1 - (pointi + 0.5)/nPoints);
            const scalar theta =
                mathematical::pi*(1 + Foam::sqrt(5.0))*(pointi + 0.5);

            hemiPts[pointi] =
                vector
                (
                    Foam::cos(theta)*Foam::sin(phi),
                    Foam::sin(theta)*Foam::sin(phi),
                    Foam::cos(phi)
                );
        }
    }
    else if (mesh_.nSolutionD() == 2)
    {
        // Point in range -1 < x < 1; y = 0; 0 < z < 1;   _\|/_

        forAll(hemiPts, pointi)
        {
            const scalar theta = mathematical::pi*(pointi+0.5)/nPoints;
            hemiPts[pointi] = vector(Foam::cos(theta), 0, Foam::sin(theta));
        }
    }

    return themiPts;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::raySearchEngine::raySearchEngine
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    mapPtr_(nullptr),
    compactToGlobal_(),
    globalNumbering_(),
    patchGroup_(dict.getOrDefault<word>("patchGroup", "viewFactorWall")),
    patchIDs_(mesh_.boundaryMesh().indices(patchGroup_)),
    patchAreas_(mesh_.boundaryMesh().nNonProcessor(), Zero),
    agglomerate_(dict.get<bool>("agglomerate")),
    agglomMeshPtr_(nullptr),
    nFace_(-1),
    nCoarseFace_(-1),
    allCf_(UPstream::nProcs()),
    allSf_(UPstream::nProcs()),
    allAgg_(UPstream::nProcs())
{
    Info<< "\nParticipating patches:" << endl;

    forAll(patchIDs_, i)
    {
        const label patchi = patchIDs_[i];
        Info<< "    " << i << ": " << mesh_.boundaryMesh()[patchi].name()
            << endl;
    }

    const word agglomName(dict.getOrDefault<word>("agglom", "finalAgglom"));

    IOobject agglomIO
    (
        agglomName,
        mesh_.facesInstance(),
        mesh_,
        IOobject::MUST_READ
    );


    if (agglomerate_)
    {
        // Sets allCf_, allSf_, allAgg_ based on coarse mesh
        // Sets nFace_, nCoarseFace_
        createAgglomeration(agglomIO);
    }
    else
    {
        // Check for presence of finalAgglom - will cause problems in later
        // calculations with viewFactor radiation model
        if (agglomIO.typeHeaderOk<labelListIOList>())
        {
            WarningInFunction
                << "Found agglomeration file: " << agglomIO.objectPath() << nl
                << "    This is inconsistent with the view factor calculation "
                << "and should be removed" << nl << endl;
        }

        // Sets allCf_, allSf_ based on fine mesh
        // Sets nFace_; nCoarseFace_ remains unset (-1)
        createGeometry();
    }

    globalNumbering_ =
        nCoarseFace_ == -1 ? globalIndex(nFace_) : globalIndex(nCoarseFace_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VF::raySearchEngine::correct
(
    labelListList& visibleFaceFaces
) const
{
    labelList rayStartFace;
    labelList rayEndFace;
    shootRays(rayStartFace, rayEndFace);

    const label nFace = nParticipatingFaces();

    // Calculate number of visible faces from each local start face
    labelList nVisibleFaceFaces(nFace, Zero);
    for (const label facei : rayStartFace)
    {
        ++nVisibleFaceFaces[facei];
    }

    check(nVisibleFaceFaces);

    createParallelAddressing(rayEndFace);

    // Set visible face-faces

    // visibleFaceFaces has:
    //    (local face, local viewed face) = compact viewed face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    visibleFaceFaces.resize_nocopy(nFace);
    forAll(nVisibleFaceFaces, facei)
    {
        visibleFaceFaces[facei].resize_nocopy(nVisibleFaceFaces[facei]);
    }

    nVisibleFaceFaces = 0;
    forAll(rayStartFace, i)
    {
        const label facei = rayStartFace[i];
        const label sloti = rayEndFace[i];
        visibleFaceFaces[facei][nVisibleFaceFaces[facei]++] = sloti;
    }
}


void Foam::VF::raySearchEngine::compactAddressing
(
    const mapDistribute& map,
    pointField& compactCf,
    vectorField& compactSf,
    List<List<vector>>& compactFineSf,
    List<List<point>>& compactFineCf,
    DynamicList<List<point>>& compactPoints,
    DynamicList<label>& compactPatchId
) const
{
    compactCf.resize_nocopy(map.constructSize());
    compactSf.resize_nocopy(map.constructSize());
    compactFineSf.resize_nocopy(map.constructSize());
    compactFineCf.resize_nocopy(map.constructSize());
    compactPoints.setCapacity(map.constructSize());
    compactPatchId.setCapacity(map.constructSize());

    // Insert my local values area and centre values
    if (agglomMeshPtr_)
    {
        SubList<vector>(compactSf, nCoarseFace_) = allSf_[Pstream::myProcNo()];
        SubList<point>(compactCf, nCoarseFace_) = allCf_[Pstream::myProcNo()];

        const auto& coarseMesh = agglomMeshPtr_();
        const auto& coarsePatches = coarseMesh.boundaryMesh();
        const auto& coarseFaces = coarseMesh.faces();
        const auto& coarsePoints = coarseMesh.points();

        const auto& finalAgglom = coarseMesh.patchFaceAgglomeration();

        // Insert my fine local values per coarse face
        label sloti = 0;
        for (const label patchi : patchIDs_)
        {
            const auto& fineCfp = mesh_.Cf().boundaryField()[patchi];
            const auto& fineSfp = mesh_.Sf().boundaryField()[patchi];
            const labelList& agglom = finalAgglom[patchi];

            if (agglom.empty()) continue;

            const label nAgglom = max(agglom) + 1;
            const labelListList coarseToFine = invertOneToMany(nAgglom, agglom);
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchi];

            const label coarseStart = coarsePatches[patchi].start();

            forAll(coarseToFine, coarsei)
            {
                compactPatchId.push_back(patchi);

                const vectorField localPoints
                (
                    coarsePoints,
                    coarseFaces[coarseStart + coarsei]
                );
                compactPoints.push_back(localPoints);

                const label coarseFacei = coarsePatchFace[coarsei];
                const labelList& fineFaces = coarseToFine[coarseFacei];

                List<point>& fineCf = compactFineCf[sloti];
                fineCf.resize_nocopy(fineFaces.size());
                fineCf = UIndirectList<point>(fineCfp, fineFaces);

                List<vector>& fineSf = compactFineSf[sloti];
                fineSf.resize_nocopy(fineFaces.size());
                fineSf = UIndirectList<vector>(fineSfp, fineFaces);

                ++sloti;
            }
        }
    }
    else
    {
        SubList<vector>(compactSf, nFace_) = allSf_[Pstream::myProcNo()];
        SubList<point>(compactCf, nFace_) = allCf_[Pstream::myProcNo()];

        const auto& patches = mesh_.boundaryMesh();
        const faceList& faces = mesh_.faces();

        label sloti = 0;

        for (const label patchi : patchIDs_)
        {
            const auto& Sfp = mesh_.Sf().boundaryField()[patchi];
            const auto& Cfp = mesh_.Cf().boundaryField()[patchi];

            const polyPatch& pp = patches[patchi];

            forAll(pp, facei)
            {
                compactPatchId.push_back(patchi);

                const auto& fpts = faces[facei + pp.start()];
                compactPoints.push_back(List<point>(mesh_.points(), fpts));

                compactFineCf[sloti] = List<point>({Cfp[facei]});
                compactFineSf[sloti] = List<vector>({Sfp[facei]});
                ++sloti;
            }
        }
    }


    // Do all swapping
    map.distribute(compactSf);
    map.distribute(compactCf);
    map.distribute(compactFineCf);
    map.distribute(compactFineSf);
    map.distribute(compactPoints);
    map.distribute(compactPatchId);
}


// ************************************************************************* //
