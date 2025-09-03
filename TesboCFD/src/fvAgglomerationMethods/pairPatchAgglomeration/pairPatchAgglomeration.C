/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "pairPatchAgglomeration.H"
#include "meshTools.H"
#include "edgeHashes.H"
#include "unitConversion.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pairPatchAgglomeration, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairPatchAgglomeration::compactLevels(const label nCreatedLevels)
{
    nFaces_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    patchLevels_.setSize(nCreatedLevels);
}


bool Foam::pairPatchAgglomeration::continueAgglomerating
(
    const label nLocal,
    const label nLocalOld,
    const label nMarkedEdges
)
{
    // Keep agglomerating
    // - if global number of faces is still changing
    // - and if local number of faces still too large (on any processor)
    //       or if global number of faces still too large

    label nGlobal = returnReduce(nLocal, sumOp<label>());
    label nGlobalOld = returnReduce(nLocalOld, sumOp<label>());
    label nGlobalMarked = returnReduce(nMarkedEdges, sumOp<label>());

    return
    (
        returnReduceOr(nLocal > nFacesInCoarsestLevel_)
     || nGlobal > nGlobalFacesInCoarsestLevel_
    )
    && (nGlobal != nGlobalOld || nGlobalMarked > 0);
}


void Foam::pairPatchAgglomeration::setLevel0EdgeWeights()
{
    const bPatch& coarsePatch = patchLevels_[0];
    const auto& coarseEdges = coarsePatch.edges();

    // Statistics on edges
    label nNonManif = 0;
    label nFeat = 0;

    for (label i = 0; i < coarsePatch.nInternalEdges(); i++)
    {
        scalar edgeLength = coarseEdges[i].mag(coarsePatch.localPoints());

        const labelList& eFaces = coarsePatch.edgeFaces()[i];

        if (eFaces.size() == 2)
        {
            scalar cosI =
                coarsePatch.faceNormals()[eFaces[0]]
              & coarsePatch.faceNormals()[eFaces[1]];

            const edge edgeCommon = edge(eFaces[0], eFaces[1]);

            if (facePairWeight_.found(edgeCommon))
            {
                facePairWeight_[edgeCommon] += edgeLength;
            }
            else
            {
                facePairWeight_.insert(edgeCommon, edgeLength);
            }

            if (cosI < Foam::cos(degToRad(featureAngle_)))
            {
                facePairWeight_[edgeCommon] = -1.0;
                nFeat++;
            }
        }
        else
        {
            forAll(eFaces, j)
            {
                for (label k = j+1; k<eFaces.size(); k++)
                {
                    facePairWeight_.insert
                    (
                        edge(eFaces[j], eFaces[k]),
                        -1.0
                    );
                }
            }
            nNonManif++;
        }
    }

    if (debug)
    {
        Pout<< "Level:" << 0
            << " nEdges:" << coarsePatch.nEdges() << " of which:" << nl
            << "    boundary:" << coarsePatch.nBoundaryEdges() << nl
            << "    non-manifold:" << nNonManif << nl
            << "    feature (angle < " << featureAngle_ << "):" << nFeat << nl
            << endl;
    }
}


void Foam::pairPatchAgglomeration::setEdgeWeights
(
    const label fineLevelIndex
)
{
    const bPatch& coarsePatch = patchLevels_[fineLevelIndex];
    const auto& coarseEdges = coarsePatch.edges();
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    edgeHashSet fineFeaturedFaces(coarsePatch.nEdges()/10);

    // Map fine faces with featured edge into coarse faces
    forAllConstIters(facePairWeight_, iter)
    {
        if (iter() == -1.0)
        {
            const edge e = iter.key();
            const edge edgeFeatured
            (
                fineToCoarse[e[0]],
                fineToCoarse[e[1]]
            );
            fineFeaturedFaces.insert(edgeFeatured);
        }
    }


    // Statistics on edges
    label nNonManif = 0;
    label nFeat = 0;

    // Clean old weights
    facePairWeight_.clear();
    facePairWeight_.reserve(coarsePatch.nEdges());

    for (label i = 0; i < coarsePatch.nInternalEdges(); i++)
    {
        scalar edgeLength = coarseEdges[i].mag(coarsePatch.localPoints());

        const labelList& eFaces = coarsePatch.edgeFaces()[i];

        if (eFaces.size() == 2)
        {
            const edge edgeCommon(eFaces[0], eFaces[1]);
            // If the fine 'pair' faces was featured edge so it is
            // the coarse 'pair'
            if (fineFeaturedFaces.found(edgeCommon))
            {
                auto w = facePairWeight_.find(edgeCommon);
                if (!w.good())
                {
                    facePairWeight_.insert(edgeCommon, -1.0);
                    nFeat++;
                }
                else if (w() != -1.0)
                {
                    // Mark as feature edge
                    w() = -1.0;
                    nFeat++;
                }
            }
            else
            {
                auto w = facePairWeight_.find(edgeCommon);
                if (w)
                {
                    if (w() != -1.0)
                    {
                        w() += edgeLength;
                    }
                }
                else
                {
                    facePairWeight_.insert(edgeCommon, edgeLength);
                }
            }
        }
        else
        {
            // Set edge as barrier by setting weight to -1
            forAll(eFaces, j)
            {
                for (label k = j+1; k<eFaces.size(); k++)
                {
                    facePairWeight_.insert
                    (
                        edge(eFaces[j], eFaces[k]),
                        -1.0
                    );
                }
            }
            nNonManif++;
        }
    }

    if (debug)
    {
        Pout<< "Level:" << fineLevelIndex
            << " nEdges:" << coarsePatch.nEdges() << " of which:" << nl
            << "    boundary:" << coarsePatch.nBoundaryEdges() << nl
            << "    non-manifold:" << nNonManif << nl
            << "    feature (angle < " << featureAngle_ << "):" << nFeat << nl
            << endl;
    }
}


bool Foam::pairPatchAgglomeration::isSingleEdgeLoop
(
    const bPatch& patch,
    const labelList& faceIDs,
    const label facei
) const
{
    // Does combining facei with faceIDs produce a valid single face?

    labelList allFaces(faceIDs.size()+1);
    SubList<label>(allFaces, faceIDs.size()) = faceIDs;
    allFaces.last() = facei;

    // Construct single face
    const indirectPrimitivePatch upp
    (
        IndirectList<face>(patch, allFaces),
        patch.points()
    );

    return (upp.edgeLoops().size() == 1);
}


Foam::label Foam::pairPatchAgglomeration::maxValidNeighbour
(
    const bool addToCluster,
    const bPatch& patch,
    const label facei,
    const labelList& fineToCoarse
    //const labelListList& coarseToFine
) const
{
    // Return index of neighbour face with max edge weight. Either looks
    // at clustered faces (addToCluster=true) or at unclustered faces

    const auto& fFaces = patch.faceFaces()[facei];

    label matchFaceNeibNo = -1;
    scalar maxFaceWeight = -0.5;    // negative but larger than -1 (= feature)

    if (addToCluster)
    {
        // Check faces to find grouped neighbour with largest face weight
        // and that forms a single edge cut
        for (const label faceNeig : fFaces)
        {
            const label coarsei = fineToCoarse[faceNeig];

            if (coarsei >= 0)
            {
                const edge edgeCommon = edge(facei, faceNeig);
                const auto& weight = facePairWeight_[edgeCommon];
                if
                (
                    (weight > maxFaceWeight)
                 //&& (isSingleEdgeLoop(patch, coarseToFine[coarsei], facei))
                )
                {
                    maxFaceWeight = weight;
                    matchFaceNeibNo = faceNeig;
                }
            }
        }
    }
    else
    {
        // Check faces to find ungrouped neighbour with largest face weight
        for (const label faceNeig : fFaces)
        {
            const label coarsei = fineToCoarse[faceNeig];

            if (coarsei < 0)    // ungrouped
            {
                const edge edgeCommon = edge(facei, faceNeig);
                const auto& weight = facePairWeight_[edgeCommon];
                if (weight > maxFaceWeight)
                {
                    maxFaceWeight = weight;
                    matchFaceNeibNo = faceNeig;
                }
            }
        }
    }

    return matchFaceNeibNo;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPatchAgglomeration::pairPatchAgglomeration
(
    const faceList& faces,
    const pointField& points,
    const dictionary& controlDict
)
:
    mergeLevels_
    (
        controlDict.getOrDefault<label>("mergeLevels", 2)
    ),
    maxLevels_(50),
    nFacesInCoarsestLevel_
    (
        controlDict.get<label>("nFacesInCoarsestLevel")
    ),
    nGlobalFacesInCoarsestLevel_(labelMax),
    //(
    //    controlDict.get<label>("nGlobalFacesInCoarsestLevel")
    //),
    featureAngle_
    (
        controlDict.getOrDefault<scalar>("featureAngle", 0)
    ),
    nFaces_(maxLevels_),
    restrictAddressing_(maxLevels_),
    restrictTopBottomAddressing_(identity(faces.size())),
    patchLevels_(maxLevels_),
    facePairWeight_(faces.size())
{
    // Set base fine patch
    patchLevels_.set(0, new bPatch(faces, points));

    // Set number of faces for the base patch
    nFaces_[0] = faces.size();

    // Set edge weights for level 0
    setLevel0EdgeWeights();
}


Foam::pairPatchAgglomeration::pairPatchAgglomeration
(
    const faceList& faces,
    const pointField& points,
    const label mergeLevels,
    const label maxLevels,
    const label nFacesInCoarsestLevel,          // local number of cells
    const label nGlobalFacesInCoarsestLevel,    // global number of cells
    const scalar featureAngle
)
:
    mergeLevels_(mergeLevels),
    maxLevels_(maxLevels),
    nFacesInCoarsestLevel_(nFacesInCoarsestLevel),
    nGlobalFacesInCoarsestLevel_(nGlobalFacesInCoarsestLevel),
    featureAngle_(featureAngle),
    nFaces_(maxLevels_),
    restrictAddressing_(maxLevels_),
    restrictTopBottomAddressing_(identity(faces.size())),
    patchLevels_(maxLevels_),
    facePairWeight_(faces.size())
{
    // Set base fine patch
    patchLevels_.set(0, new bPatch(faces, points));

    // Set number of faces for the base patch
    nFaces_[0] = faces.size();

    // Set edge weights for level 0
    setLevel0EdgeWeights();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairPatchAgglomeration::~pairPatchAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pairPatchAgglomeration::bPatch&
Foam::pairPatchAgglomeration::patchLevel
(
    const label i
) const
{
    return patchLevels_[i];
}


void Foam::pairPatchAgglomeration::mapBaseToTopAgglom
(
    const label fineLevelIndex
)
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];
    forAll(restrictTopBottomAddressing_, i)
    {
        restrictTopBottomAddressing_[i] =
            fineToCoarse[restrictTopBottomAddressing_[i]];
    }
}


bool Foam::pairPatchAgglomeration::agglomeratePatch
(
    const bPatch& patch,
    const labelList& fineToCoarse,
    const label fineLevelIndex,
    label& nMarkedEdges
)
{
    if (min(fineToCoarse) == -1)
    {
        FatalErrorInFunction
            << "min(fineToCoarse) == -1" << exit(FatalError);
    }

    if (fineToCoarse.size() == 0)
    {
        return false;
    }

    if (fineToCoarse.size() != patch.size())
    {
        FatalErrorInFunction
            << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << fineToCoarse.size()
            << " nEqns: " << patch.size()
            << abort(FatalError);
    }

    const label nCoarseI =  max(fineToCoarse) + 1;
    List<face> patchFaces(nCoarseI);

    // Patch faces per agglomeration
    labelListList coarseToFine(invertOneToMany(nCoarseI, fineToCoarse));

    // Additional feature edges created
    nMarkedEdges = 0;


    for (label coarseI = 0; coarseI < nCoarseI; coarseI++)
    {
        const labelList& fineFaces = coarseToFine[coarseI];

        // Construct single face
        indirectPrimitivePatch upp
        (
            IndirectList<face>(patch, fineFaces),
            patch.points()
        );

        if (upp.edgeLoops().size() != 1)
        {
            // More than one outside loop. Possible because pair-wise
            // agglomeration has e.g. walked a path leaving a hole in the
            // middle of a coarse face.
            //  - mark as feature edges to avoid locally further agglomeration
            //  - but further agglomeration might e.g. 'fill in the hole'
            //  - or immediately leave all agglomeration
            // Currently we have no choice but to leave agglomeration since
            // we cannot store a face-with-hole.


            //{
            //    OBJstream os
            //    (
            //        "error_agglomeration_"+Foam::name(fineLevelIndex)+".obj"
            //    );
            //    Pout<< "Writing error patch at level:" << fineLevelIndex
            //        << " to:" << os.name() << endl;
            //    os.write(upp.localFaces(), upp.localPoints(), true);
            //}

            if (fineFaces.size() >= 2)
            {
                forAll(fineFaces, j)
                {
                    for (label k = j+1; k<fineFaces.size(); k++)
                    {
                        const edge e(fineFaces[j], fineFaces[k]);

                        auto w = facePairWeight_.find(e);
                        if (!w.good())
                        {
                            facePairWeight_.insert(e, -1.0);
                            nMarkedEdges++;
                        }
                        else if (w() != -1.0)
                        {
                            // Mark as feature edge
                            w() = -1.0;
                            nMarkedEdges++;
                        }
                    }
                }
            }

            return false;
        }

        // In-place override face
        patchFaces[coarseI] = face
        (
            renumber
            (
                upp.meshPoints(),
                upp.edgeLoops()[0]
            )
        );
    }

    patchLevels_.set
    (
        fineLevelIndex,
        new bPatch
        (
            SubList<face>(patchFaces, nCoarseI, 0),
            patch.points()
        )
    );

    return true;
}


void Foam::pairPatchAgglomeration::agglomerate()
{
    label nPairLevels = 0;
    label nCreatedLevels = 1; // 0 level is the base patch

    label nCoarseFaces = 0;
    label nCoarseFacesOld = 0;
    label nMarkedEdges = 0;

    while (nCreatedLevels < maxLevels_)
    {
        const bPatch& patch = patchLevels_[nCreatedLevels - 1];

        // Agglomerate locally
        tmp<labelField> tfinalAgglom;

        bool createdLevel = false;
        while (!createdLevel)
        {
            // Agglomerate locally using edge weights
            // - calculates nCoarseFaces; returns fine to coarse addressing
            tfinalAgglom = agglomerateOneLevel(nCoarseFaces, patch);

            if (nCoarseFaces == 0)
            {
                break;
            }
            else
            {
                // Attempt to create coarse face addressing
                // - returns true if successful; otherwise resets edge weights
                //   and tries again...
                createdLevel = agglomeratePatch
                (
                    patch,
                    tfinalAgglom,
                    nCreatedLevels,
                    nMarkedEdges
                );
            }
        }

        if (createdLevel)
        {
            if (debug)
            {
                const auto& agglomPatch = patchLevels_[nCreatedLevels];
                OBJstream os("agglomPatch"+Foam::name(nCreatedLevels)+".obj");
                Pout<< "Writing new patch at level:" << nCreatedLevels
                    << " to:" << os.name() << endl;
                os.write(agglomPatch, agglomPatch.points(), true);
            }

            restrictAddressing_.set(nCreatedLevels, tfinalAgglom);

            mapBaseToTopAgglom(nCreatedLevels);

            setEdgeWeights(nCreatedLevels);

            if (nPairLevels % mergeLevels_)
            {
                combineLevels(nCreatedLevels);
            }
            else
            {
                nCreatedLevels++;
            }

            nPairLevels++;

            nFaces_[nCreatedLevels] = nCoarseFaces;
        }

        // Check to see if we need to continue agglomerating
        // - Note: performs parallel reductions
        if (!continueAgglomerating(nCoarseFaces, nCoarseFacesOld, nMarkedEdges))
        {
            break;
        }

        nCoarseFacesOld = nCoarseFaces;
    }

    compactLevels(nCreatedLevels);
}


Foam::tmp<Foam::labelField> Foam::pairPatchAgglomeration::agglomerateOneLevel
(
    label& nCoarseFaces,
    const bPatch& patch
)
{
    const label nFineFaces = patch.size();

    auto tcoarseCellMap = tmp<labelField>::New(nFineFaces, -1);
    auto& coarseCellMap = tcoarseCellMap.ref();

    const labelListList& faceFaces = patch.faceFaces();

    //labelListList coarseToFine(nFineFaces);
    nCoarseFaces = 0;

    forAll(faceFaces, facei)
    {
        if (coarseCellMap[facei] < 0)
        {
            const label matchFaceNeibNo = maxValidNeighbour
            (
                false,          // ungrouped neighbours only
                patch,
                facei,
                coarseCellMap
                //coarseToFine
            );

            if (matchFaceNeibNo >= 0)
            {
                // Make a new group
                coarseCellMap[facei] = nCoarseFaces;
                coarseCellMap[matchFaceNeibNo] = nCoarseFaces;
                //coarseToFine[nCoarseFaces] =
                //    labelList({facei, matchFaceNeibNo});
                nCoarseFaces++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                const label clusterMatchFaceNo = maxValidNeighbour
                (
                    true,          // grouped neighbours only
                    patch,
                    facei,
                    coarseCellMap
                    //coarseToFine
                );

                if (clusterMatchFaceNo >= 0)
                {
                    // Add the cell to the best cluster
                    const label coarsei = coarseCellMap[clusterMatchFaceNo];
                    coarseCellMap[facei] = coarsei;
                    //coarseToFine[coarsei].append(facei);
                }
                else
                {
                    // If not create single-cell "clusters" for each
                    coarseCellMap[facei] = nCoarseFaces;
                    //coarseToFine[nCoarseFaces] = labelList({facei});
                    nCoarseFaces++;
                }
            }
        }
    }

    //coarseToFine.setSize(nCoarseFaces);

    // Check that all faces are part of clusters,
    for (label facei=0; facei<nFineFaces; facei++)
    {
        if (coarseCellMap[facei] < 0)
        {
            FatalErrorInFunction
                << " face " << facei
                << " is not part of a cluster"
                << exit(FatalError);
        }
    }

    return tcoarseCellMap;
}


void Foam::pairPatchAgglomeration::combineLevels(const label curLevel)
{
    label prevLevel = curLevel - 1;

    // Set the previous level nCells to the current
    nFaces_[prevLevel] = nFaces_[curLevel];

    // Map the restrictAddressing from the coarser level into the previous
    // finer level

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    // Delete the restrictAddressing for the coarser level
    restrictAddressing_.set(curLevel, nullptr);

    patchLevels_.set(prevLevel, patchLevels_.set(curLevel, nullptr));
}


// ************************************************************************* //
