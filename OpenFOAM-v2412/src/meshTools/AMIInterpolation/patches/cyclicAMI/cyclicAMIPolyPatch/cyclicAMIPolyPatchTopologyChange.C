/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "cyclicAMIPolyPatch.H"
#include "SubField.H"
#include "vectorList.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::restoreScaledGeometry()
{
    DebugInFunction << endl;

    // Note: only used for topology update (createAMIFaces_ flag)
    if (!createAMIFaces_)
    {
        FatalErrorInFunction
            << "Attempted to perform topology update when createAMIFaces_ "
            << "flag is set to false"
            << abort(FatalError);
    }

    if (boundaryMesh().mesh().hasCellVolumes())
    {
        WarningInFunction
            << "Mesh already has volumes set!"
            << endl;
    }

    vectorField::subField faceAreas = this->faceAreas();
    vectorField::subField faceCentres = this->faceCentres();

    DebugInfo
        << "Patch:" << name() << " before: sum(mag(faceAreas)):"
        << gSum(mag(faceAreas)) << nl
        << "Patch:" << name() << " before: sum(mag(faceAreas0)):"
        << gSum(mag(faceAreas0_)) << endl;

    faceAreas = faceAreas0_;
    if (moveFaceCentres_)
    {
        DebugInfo << "Moving face centres" << endl;
        faceCentres = faceCentres0_;
    }

    faceAreas0_.clear();
    faceCentres0_.clear();

    DebugInfo
        << "Patch:" << name() << " after: sum(mag(faceAreas)):"
        << gSum(mag(faceAreas)) << nl
        << "Patch:" << name() << " after: sum(mag(faceAreas0)):"
        << gSum(mag(faceAreas0_)) << endl;
}


bool Foam::cyclicAMIPolyPatch::removeAMIFaces(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    // Note: only used for topology update (createAMIFaces_ flag)
    if (!createAMIFaces_)
    {
        FatalErrorInFunction
            << "Attempted to perform topology update when createAMIFaces_ "
            << "flag is set to false"
            << abort(FatalError);
    }

    if (!owner())
    {
        return false;
    }

    bool changeRequired = false;

    // Remove any faces that we inserted to create the 1-to-1 match...

    const cyclicAMIPolyPatch& nbr = neighbPatch();

    const label newSrcFaceStart = srcFaceIDs_.size();

    if (newSrcFaceStart != 0)
    {
        for (label facei = newSrcFaceStart; facei < size(); ++facei)
        {
            changeRequired = true;
            label meshFacei = start() + facei;
            topoChange.removeFace(meshFacei, -1);
        }
    }

    const label newTgtFaceStart = tgtFaceIDs_.size();

    if (newTgtFaceStart != 0)
    {
        for (label facei = newTgtFaceStart; facei < nbr.size(); ++facei)
        {
            changeRequired = true;
            label meshFacei = nbr.start() + facei;
            topoChange.removeFace(meshFacei, -1);
        }
    }

    srcFaceIDs_.clear();
    tgtFaceIDs_.clear();

    return changeRequired;
}


bool Foam::cyclicAMIPolyPatch::addAMIFaces(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    // Note: only used for topology update (createAMIFaces_ flag = true)
    if (!createAMIFaces_)
    {
        FatalErrorInFunction
            << "Attempted to perform topology update when createAMIFaces_ "
            << "flag is set to false"
            << abort(FatalError);
    }

    bool changedFaces = false;
    const cyclicAMIPolyPatch& nbr = neighbPatch();

    polyMesh& mesh = const_cast<polyMesh&>(boundaryMesh().mesh());
    const faceZoneMesh& faceZones = mesh.faceZones();

    // First face address and weight are used to manipulate the
    // original face - all other addresses and weights are used to
    // create additional faces
    const labelListList& srcToTgtAddr = AMI().srcAddress();
    const labelListList& tgtToSrcAddr = AMI().tgtAddress();

    const label nSrcFace = srcToTgtAddr.size();
    const label nTgtFace = tgtToSrcAddr.size();

    srcFaceIDs_.setSize(nSrcFace);
    tgtFaceIDs_.setSize(nTgtFace);

    label nNewSrcFaces = 0;
    forAll(srcToTgtAddr, srcFacei)
    {
        const labelList& tgtAddr = srcToTgtAddr[srcFacei];

        // No tgt faces linked to srcFacei (ACMI)
        if (tgtAddr.empty()) continue;

        srcFaceIDs_[srcFacei].setSize(tgtAddr.size());
        srcFaceIDs_[srcFacei][0] = srcFacei;

        const label meshFacei = start() + srcFacei;
        for (label addri = 1; addri < tgtAddr.size(); ++addri)
        {
            changedFaces = true;

            // Note: new faces reuse originating face points
            // - but areas are scaled by the weights (later)

            // New source face for each target face address
            srcFaceIDs_[srcFacei][addri] = nNewSrcFaces + nSrcFace;
            ++nNewSrcFaces;
            (void)topoChange.addFace
            (
                mesh.faces()[meshFacei],        // modified face
                mesh.faceOwner()[meshFacei],    // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                meshFacei,                      // master face
                false,                          // face flip
                index(),                        // patch for face
                faceZones.whichZone(meshFacei), // zone for original face
                false                           // face flip in zone
            );
        }
    }

    label nNewTgtFaces = 0;
    forAll(tgtToSrcAddr, tgtFacei)
    {
        const labelList& srcAddr = tgtToSrcAddr[tgtFacei];

        // No src faces linked to tgtFacei (ACMI)
        if (srcAddr.empty()) continue;

        tgtFaceIDs_[tgtFacei].setSize(srcAddr.size());
        tgtFaceIDs_[tgtFacei][0] = tgtFacei;

        const label meshFacei = nbr.start() + tgtFacei;
        for (label addri = 1; addri < srcAddr.size(); ++addri)
        {
            changedFaces = true;

            // Note: new faces reuse originating face points
            // - but areas are scaled by the weights (later)

            // New target face for each source face address
            tgtFaceIDs_[tgtFacei][addri] = nNewTgtFaces + nTgtFace;
            ++nNewTgtFaces;

            (void)topoChange.addFace
            (
                mesh.faces()[meshFacei],        // modified face
                mesh.faceOwner()[meshFacei],    // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                meshFacei,                      // master face
                false,                          // face flip
                nbr.index(),                    // patch for face
                faceZones.whichZone(meshFacei), // zone for original face
                false                           // face flip in zone
            );
        }
    }

    Info<< "AMI: Patch " << name() << " additional faces: "
        << returnReduce(nNewSrcFaces, sumOp<label>()) << nl
        << "AMI: Patch " << nbr.name() << " additional faces: "
        << returnReduce(nNewTgtFaces, sumOp<label>())
        << endl;

    if (debug)
    {
        Pout<< "New faces - " << name() << ": " << nNewSrcFaces
            << " "  << nbr.name() << ": " << nNewTgtFaces << endl;
    }

    return returnReduceOr(changedFaces);
}


void Foam::cyclicAMIPolyPatch::setAMIFaces()
{
    // Create new mesh faces so that there is a 1-to-1 correspondence
    // between faces on each side of the AMI

    // Note: only used for topology update (createAMIFaces_ flag)
    if (!createAMIFaces_)
    {
        FatalErrorInFunction
            << "Attempted to perform topology update when createAMIFaces_ "
            << "flag is set to false"
            << abort(FatalError);
    }


    DebugInFunction << endl;

    if (!owner())
    {
        return;
    }

    const cyclicAMIPolyPatch& nbr = neighbPatch();

    vectorField& nbrFaceAreas0 = nbr.faceAreas0();
    vectorField& nbrFaceCentres0 = nbr.faceCentres0();

    // Scale the new face areas and set the centroids
    // Note:
    // - storing local copies so that they can be re-applied after the call to
    //   movePoints that will reset any changes to the areas and centroids
    //
    // - For AMI, src and tgt patches should be the same
    // - For ACMI they are likely to be different!
    faceAreas0_ = faceAreas();
    faceCentres0_ = faceCentres();
    nbrFaceAreas0 = nbr.faceAreas();
    nbrFaceCentres0 = nbr.faceCentres();

    // Original AMI info (based on the mesh state when the AMI was evaluated)
    const labelListList& srcToTgtAddr0 = AMIPtr_->srcAddress();
    const labelListList& tgtToSrcAddr0 = AMIPtr_->tgtAddress();
    const pointListList& srcCtr0 = AMIPtr_->srcCentroids();
    const scalarListList& srcToTgtWght0 = AMIPtr_->srcWeights();
    const label singlePatchProc = AMIPtr_->singlePatchProc();

    // New addressing on new mesh (extended by polyTopoChange)
    labelListList srcToTgtAddr1(size(), labelList());
    labelListList tgtToSrcAddr1(nbr.size(), labelList());

    // Need to calc new parallel maps (mesh has changed since AMI was computed)
    autoPtr<mapDistribute> srcToTgtMap1;
    autoPtr<mapDistribute> tgtToSrcMap1;

    if (AMIPtr_->distributed())
    {
       // Parallel running

        // Global index based on old patch sizes (when AMI was computed)
        globalIndex globalSrcFaces0(srcToTgtAddr0.size(), AMIPtr_().comm());
        globalIndex globalTgtFaces0(tgtToSrcAddr0.size(), AMIPtr_().comm());

        // Global index based on new patch sizes
        globalIndex globalSrcFaces1(size(), AMIPtr_().comm());
        globalIndex globalTgtFaces1(nbr.size(), AMIPtr_().comm());


        // Gather source side info
        // =======================

        // Note: using new global index for addressing, and distributed using
        // the old AMI map
        labelListList newTgtGlobalFaces(tgtFaceIDs_);
        forAll(newTgtGlobalFaces, tgtFacei)
        {
            globalTgtFaces1.inplaceToGlobal(newTgtGlobalFaces[tgtFacei]);
        }
        AMIPtr_->tgtMap().distribute(newTgtGlobalFaces);

        // Now have new tgt face indices for each src face

        labelList globalSrcFaceIDs(identity(srcToTgtAddr0.size()));
        globalSrcFaces0.inplaceToGlobal(globalSrcFaceIDs);
        AMIPtr_->srcMap().distribute(globalSrcFaceIDs);
        // globalSrcFaceIDs now has remote data for each srcFacei0 known to the
        // tgt patch

        List<List<point>> globalSrcCtrs0(srcCtr0);
        AMIPtr_->srcMap().distribute(globalSrcCtrs0);

        labelList globalTgtFaceIDs(identity(tgtToSrcAddr0.size()));
        globalTgtFaces0.inplaceToGlobal(globalTgtFaceIDs);
        AMIPtr_->tgtMap().distribute(globalTgtFaceIDs);
        // globalTgtFaceIDs now has remote data for each tgtFacei0 known to the
        // src patch

        // For debug - send tgt face centres and compare against mapped src
        // face centres
        //List<List<point>> globalTgtCtrs0(tgtCtr0);
        //AMIPtr_->tgtMap().distribute(globalTgtCtrs0);

        labelListList globalTgtToSrcAddr(tgtToSrcAddr0);
        forAll(tgtToSrcAddr0, tgtFacei0)
        {
            forAll(tgtToSrcAddr0[tgtFacei0], addri)
            {
                const label globalSrcFacei =
                    globalSrcFaceIDs[tgtToSrcAddr0[tgtFacei0][addri]];
                globalTgtToSrcAddr[tgtFacei0][addri] = globalSrcFacei;
            }
        }
        AMIPtr_->tgtMap().distribute(globalTgtToSrcAddr);

        labelListList globalSrcToTgtAddr(srcToTgtAddr0);
        forAll(srcToTgtAddr0, srcFacei0)
        {
            forAll(srcToTgtAddr0[srcFacei0], addri)
            {
                const label globalTgtFacei =
                    globalTgtFaceIDs[srcToTgtAddr0[srcFacei0][addri]];
                globalSrcToTgtAddr[srcFacei0][addri] = globalTgtFacei;
            }
        }
        AMIPtr_->srcMap().distribute(globalSrcToTgtAddr);

        label nError = 0;
        forAll(srcToTgtAddr0, srcFacei0)
        {
            const labelList& newSrcFaces = srcFaceIDs_[srcFacei0];

            forAll(newSrcFaces, i)
            {
                const label srcFacei1 = newSrcFaces[i];

                // What index did srcFacei0 appear in tgtToSrc0 list?
                // - if first index, all ok
                // - else tgt face has been moved to according to tgtFaceIDs_
                const label tgtFacei0 = srcToTgtAddr0[srcFacei0][i];
                const label addri =
                    globalTgtToSrcAddr[tgtFacei0].find
                    (
                        globalSrcFaceIDs[srcFacei0]
                    );

                if (addri == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find global source face "
                            << globalSrcFaceIDs[srcFacei0]
                            << " in globalTgtToSrcAddr[" << tgtFacei0 << "]: "
                            << globalTgtToSrcAddr[tgtFacei0]
                            << endl;
                    }
                }

                const label tgtFacei1 = newTgtGlobalFaces[tgtFacei0][addri];

                // Sanity check to see that we've picked the correct face
                // point tgtCtr0(globalTgtCtrs0[tgtFacei0][addri]);
                // Pout<< "srcCtr:" << srcCtr0[srcFacei0][i]
                //     << " tgtCtr:" << tgtCtr0 << endl;

                srcToTgtAddr1[srcFacei1] = labelList(1, tgtFacei1);
                faceAreas0_[srcFacei1] *= srcToTgtWght0[srcFacei0][i];
                faceCentres0_[srcFacei1] = srcCtr0[srcFacei0][i];
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError << " global source faces"
                << abort(FatalError);
        }


        // Gather Target side info
        // =======================

        labelListList newSrcGlobalFaces(srcFaceIDs_);
        forAll(newSrcGlobalFaces, srcFacei)
        {
            globalSrcFaces1.inplaceToGlobal(newSrcGlobalFaces[srcFacei]);
        }

        AMIPtr_->srcMap().distribute(newSrcGlobalFaces);

        // Now have new src face indices for each tgt face
        forAll(tgtToSrcAddr0, tgtFacei0)
        {
            const labelList& newTgtFaces = tgtFaceIDs_[tgtFacei0];
            forAll(newTgtFaces, i)
            {
                const label srcFacei0 = tgtToSrcAddr0[tgtFacei0][i];

                const label addri =
                    globalSrcToTgtAddr[srcFacei0].find
                    (
                        globalTgtFaceIDs[tgtFacei0]
                    );

                if (addri == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find global target face "
                            << globalTgtFaceIDs[tgtFacei0]
                            << " in globalSrcToTgtAddr[" << srcFacei0 << "]: "
                            << globalSrcToTgtAddr[srcFacei0]
                            << endl;
                    }
                }

                const label srcFacei1 = newSrcGlobalFaces[srcFacei0][addri];

                // Sanity check to see that we've picked the correct face
                point srcCtr0(globalSrcCtrs0[srcFacei0][addri]);
                reverseTransformPosition(srcCtr0, srcFacei0);

                const label tgtFacei1 = newTgtFaces[i];
                tgtToSrcAddr1[tgtFacei1] = labelList(1, srcFacei1);
                nbrFaceCentres0[tgtFacei1] = srcCtr0;
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError << " global target faces"
                << abort(FatalError);
        }

        // Update the maps
        {
            List<Map<label>> cMap;
            srcToTgtMap1.reset
            (
                new mapDistribute(globalSrcFaces1, tgtToSrcAddr1, cMap)
            );
        }
        {
            List<Map<label>> cMap;
            tgtToSrcMap1.reset
            (
                new mapDistribute(globalTgtFaces1, srcToTgtAddr1, cMap)
            );
        }

        // Reset tgt patch areas using the new map
        vectorList newSrcGlobalFaceAreas(faceAreas0_);

        srcToTgtMap1->distribute(newSrcGlobalFaceAreas);
        forAll(nbrFaceAreas0, tgtFacei)
        {
            if (!tgtToSrcAddr1[tgtFacei].empty())
            {
                const label srcFacei = tgtToSrcAddr1[tgtFacei][0];
                nbrFaceAreas0[tgtFacei] = -newSrcGlobalFaceAreas[srcFacei];
            }
        }
    }
    else
    {
        label nError = 0;
        forAll(srcToTgtAddr0, srcFacei0)
        {
            const labelList& srcFaceTgtAddr0 = srcToTgtAddr0[srcFacei0];
            const scalarList& srcFaceTgtWght0 = srcToTgtWght0[srcFacei0];
            const pointList& srcFaceTgtCtr0 = srcCtr0[srcFacei0];
            forAll(srcFaceTgtAddr0, addri)
            {
                const label srcFacei1 = srcFaceIDs_[srcFacei0][addri];

                // Find which slot srcFacei0 appears in tgt->src addressing
                const label tgtFacei0 = srcFaceTgtAddr0[addri];
                const label tgtAddri0 =
                    tgtToSrcAddr0[tgtFacei0].find(srcFacei0);

                if (tgtAddri0 == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find source face " << srcFacei0
                            << " in tgtToSrcAddr0[" << tgtFacei0 << "]: "
                            << tgtToSrcAddr0[tgtFacei0]
                            << endl;
                    }
                }

                const label tgtFacei1 = tgtFaceIDs_[tgtFacei0][tgtAddri0];

                faceAreas0_[srcFacei1] *= srcFaceTgtWght0[addri];
                nbrFaceAreas0[tgtFacei1] = -faceAreas0_[srcFacei1];

                point pt(srcFaceTgtCtr0[addri]);
                faceCentres0_[srcFacei1] = pt;
                reverseTransformPosition(pt, srcFacei0);
                nbrFaceCentres0[tgtFacei1] = pt;

                // SANITY CHECK
                // Info<< "srcPt:" << srcFaceCentres[srcFacei1]
                //     << " tgtPt:" << tgtFaceCentres[tgtFacei1] << endl;

                srcToTgtAddr1[srcFacei1] = labelList(1, tgtFacei1);
                tgtToSrcAddr1[tgtFacei1] = labelList(1, srcFacei1);
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError
                << " source faces in tgtToSrcAddr0"
                << abort(FatalError);
        }
    }

    scalarListList newSrcToTgtWeights(srcToTgtAddr1.size());
    forAll(srcToTgtAddr1, facei)
    {
        if (srcToTgtAddr1[facei].size())
        {
            newSrcToTgtWeights[facei] = scalarList(1, scalar(1));
        }
        else
        {
            // No connection - effect of face removed by setting area to a
            // a small value
            faceAreas0_[facei] *= tolerance_;
        }
    }

    scalarListList newTgtToSrcWeights(tgtToSrcAddr1.size());
    forAll(tgtToSrcAddr1, facei)
    {
        if (tgtToSrcAddr1[facei].size())
        {
            newTgtToSrcWeights[facei] = scalarList(1, scalar(1));
        }
        else
        {
            // No connection - effect of face removed by setting area to a
            // a small value
            nbrFaceAreas0[facei] *= tolerance_;
        }
    }

    // Reset the AMI addressing and weights to reflect the new 1-to-1
    // correspondence
    AMIPtr_->reset
    (
        std::move(srcToTgtMap1),
        std::move(tgtToSrcMap1),
        std::move(srcToTgtAddr1),
        std::move(newSrcToTgtWeights),
        std::move(tgtToSrcAddr1),
        std::move(newTgtToSrcWeights),
        singlePatchProc
    );

    // Need to set areas, e.g. for agglomeration to (re-)normalisation weights
    AMIPtr_->srcMagSf() = mag(faceAreas0_);
    AMIPtr_->tgtMagSf() = mag(nbrFaceAreas0);

    if (debug)
    {
        Pout<< "cyclicAMIPolyPatch : " << name()
            << " constructed AMI with " << nl
            << "    " << "srcAddress:" << AMIPtr_().srcAddress().size()
            << nl
            << "    " << "tgAddress :" << AMIPtr_().tgtAddress().size()
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicAMIPolyPatch::changeTopology() const
{
    DebugInFunction << endl;

    createAMIFaces_ = true;

    return true;
}


bool Foam::cyclicAMIPolyPatch::setTopology(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    if (createAMIFaces_ && owner())
    {
        // Calculate the AMI using the new points
        // Note: mesh still has old points
        resetAMI(topoChange.points());

        removeAMIFaces(topoChange);

        addAMIFaces(topoChange);

        return true;
    }

    return false;
}


// ************************************************************************* //
