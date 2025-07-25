/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019,2023 OpenCFD Ltd.
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

#include "AMIInterpolation.H"
#include "cyclicACMIGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicACMIGAMGInterface,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicACMIGAMGInterface,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterface::cyclicACMIGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces
    ),
    neighbPatchID_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).neighbPatchID()
    ),
    owner_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).owner()
    ),
    forwardT_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).forwardT()
    ),
    reverseT_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).reverseT()
    )
{
    const auto& fineCyclicACMIInterface =
        refCast<const cyclicACMILduInterface>(fineInterface);

    // Construct face agglomeration from cell agglomeration
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        for (const label curMaster : localRestrictAddressing)
        {
            const auto iter = masterToCoarseFace.cfind(curMaster);

            if (iter.good())
            {
                // Already have coarse face
                dynFaceRestrictAddressing.append(iter.val());
            }
            else
            {
                // New coarse face
                const label coarseI = dynFaceCells.size();
                dynFaceRestrictAddressing.append(coarseI);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseI);
            }
        }

        faceCells_.transfer(dynFaceCells);
        faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
    }


    // On the owner side construct the AMI

    if (fineCyclicACMIInterface.owner())
    {
        // Construct the neighbour side agglomeration (as the neighbour would
        // do it so it the exact loop above using neighbourRestrictAddressing
        // instead of localRestrictAddressing)

        labelList nbrFaceRestrictAddressing;
        {
            // From face to coarse face
            DynamicList<label> dynNbrFaceRestrictAddressing
            (
                neighbourRestrictAddressing.size()
            );

            Map<label> masterToCoarseFace(neighbourRestrictAddressing.size());

            for (const label curMaster : neighbourRestrictAddressing)
            {
                const auto iter = masterToCoarseFace.cfind(curMaster);

                if (iter.good())
                {
                    // Already have coarse face
                    dynNbrFaceRestrictAddressing.append(iter.val());
                }
                else
                {
                    // New coarse face
                    const label coarseI = masterToCoarseFace.size();
                    dynNbrFaceRestrictAddressing.append(coarseI);
                    masterToCoarseFace.insert(curMaster, coarseI);
                }
            }

            nbrFaceRestrictAddressing.transfer(dynNbrFaceRestrictAddressing);
        }


        amiPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                fineCyclicACMIInterface.AMI(),
                faceRestrictAddressing_,
                nbrFaceRestrictAddressing
            )
        );
    }
}


Foam::cyclicACMIGAMGInterface::cyclicACMIGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    GAMGInterface(index, coarseInterfaces, is),
    neighbPatchID_(readLabel(is)),
    owner_(readBool(is)),
    forwardT_(is),
    reverseT_(is)
{
    const bool hasAMI(readBool(is));

    if (hasAMI)
    {
        amiPtr_.reset(new AMIPatchToPatchInterpolation(is));
    }
}


Foam::cyclicACMIGAMGInterface::cyclicACMIGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelList& interfaceMap,
    const labelUList& faceCells,
    const labelUList& faceRestrictAddresssing,
    const labelUList& faceOffsets,
    const lduInterfacePtrsList& allInterfaces,
    const label coarseComm,
    const label myProcNo,
    const labelList& procAgglomMap
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces,
        faceCells,
        faceRestrictAddresssing
    ),
    neighbPatchID_
    (
        interfaceMap.find
        (
            refCast
            <
                const cyclicACMILduInterface
            >(fineInterface).neighbPatchID()
        )
    ),
    owner_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).owner()
    ),
    forwardT_
    (
        refCast<const cyclicACMILduInterface>(fineInterface).forwardT()
    ),
    reverseT_
    (
       refCast<const cyclicACMILduInterface>(fineInterface).reverseT()
    )
{
    const auto& fineCyclicACMIInterface =
        refCast<const cyclicACMIGAMGInterface>(fineInterface);

    if (fineCyclicACMIInterface.amiPtr_)
    {
        const auto& AMI = const_cast<AMIPatchToPatchInterpolation&>
        (
            fineCyclicACMIInterface.AMI()
        );

        label singlePatchProc = AMI.singlePatchProc();


        // Get some sizes
        label nSrc = 0;
        label nTgt = 0;
        bool hasSrcMagSf = false;
        bool hasSrcCentroids = false;
        bool hasTgtMagSf = false;

        forAll(allInterfaces, inti)
        {
            if (allInterfaces.set(inti))
            {
                const auto& intf =
                    refCast<const cyclicACMIGAMGInterface>(allInterfaces[inti]);
                const auto& AMI = intf.AMI();
                nSrc += AMI.srcAddress().size();
                nTgt += AMI.tgtAddress().size();

                if (AMI.srcMagSf().size())
                {
                    hasSrcMagSf = true;
                    if (AMI.srcMagSf().size() != AMI.srcAddress().size())
                    {
                        FatalErrorInFunction
                            << "srcMagSf size:" << AMI.srcMagSf().size()
                            << "srcAddress size:" << AMI.srcAddress().size()
                            << exit(FatalError);
                    }
                }
                if (AMI.srcCentroids().size())
                {
                    hasSrcCentroids = true;
                    if (AMI.srcCentroids().size() != AMI.srcAddress().size())
                    {
                        FatalErrorInFunction
                            << "srcCentroids size:" << AMI.srcCentroids().size()
                            << "srcAddress size:" << AMI.srcAddress().size()
                            << exit(FatalError);
                    }
                }
                if (AMI.tgtMagSf().size())
                {
                    hasTgtMagSf = true;
                    if (AMI.tgtMagSf().size() != AMI.tgtAddress().size())
                    {
                        FatalErrorInFunction
                            << "tgtMagSf size:" << AMI.tgtMagSf().size()
                            << "tgtAddress size:" << AMI.tgtAddress().size()
                            << exit(FatalError);
                    }
                }
            }
        }


        labelListList srcAddress;
        scalarListList srcWeights;
        scalarList srcMagSf;
        // Needed?
        pointListList srcCentroids;

        labelListList tgtAddress;
        scalarListList tgtWeights;
        scalarList tgtMagSf;


        // Map to send src side data to tgt side
        autoPtr<mapDistribute> srcToTgtMap;

        // Map to send tgt side data to src side
        autoPtr<mapDistribute> tgtToSrcMap;

        if (AMI.distributed())
        {
            // Create combined maps
            UPtrList<const mapDistribute> srcMaps(allInterfaces.size());
            UPtrList<const mapDistribute> tgtMaps(allInterfaces.size());
            forAll(allInterfaces, inti)
            {
                if (allInterfaces.set(inti))
                {
                    const auto& intf = refCast<const cyclicACMIGAMGInterface>
                    (
                        allInterfaces[inti]
                    );
                    const auto& AMI = intf.AMI();
                    srcMaps.set(inti, &AMI.srcMap());
                    tgtMaps.set(inti, &AMI.tgtMap());
                }
            }


            // Find ranks that agglomerate together
            const label myAgglom =
                procAgglomMap[UPstream::myProcNo(AMI.comm())];

            // Invert procAgglomMap
            const labelListList newToOldRanks
            (
                invertOneToMany
                (
                    UPstream::nProcs(coarseComm),
                    procAgglomMap
                )
            );
            const labelList& localRanks = newToOldRanks[myAgglom];


            // Offsets for slots into results of srcToTgtMap
            labelList srcStartOfLocal;
            List<Map<label>> srcCompactMaps;

            srcToTgtMap.reset
            (
                new mapDistribute
                (
                    srcMaps,
                    localRanks,     // per src map which rank it is from
                    coarseComm,
                    newToOldRanks,  // destination rank to source ranks
                    srcStartOfLocal,
                    srcCompactMaps
                )
            );

            // Assemble tgtAddress
            tgtAddress.setSize(nTgt);
            if (tgtAddress.size())
            {
                label alli = 0;
                forAll(allInterfaces, inti)
                {
                    if (allInterfaces.set(inti))
                    {
                        const auto& intf =
                            refCast<const cyclicACMIGAMGInterface>
                            (
                                allInterfaces[inti]
                            );
                        const auto& AMI = intf.AMI();
                        const auto& tgtSlots = AMI.tgtAddress();
                        const label localSize =
                            srcStartOfLocal[inti+1]
                          - srcStartOfLocal[inti];

                        forAll(tgtSlots, tgti)
                        {
                            // Append old slots: copy old values and adapt
                            auto& newSlots = tgtAddress[alli++];
                            newSlots = tgtSlots[tgti];

                            // Renumber to new indices
                            mapDistributeBase::renumberMap
                            (
                                newSlots,
                                localSize,
                                srcStartOfLocal[inti],
                                srcCompactMaps[inti],
                                AMI.srcMap().constructHasFlip() //hasFlip
                            );

                            for (const label slot : newSlots)
                            {
                                if
                                (
                                    slot < 0
                                 || slot >= srcToTgtMap().constructSize()
                                )
                                {
                                    FatalErrorInFunction << " newSlots:"
                                        << newSlots << exit(FatalError);
                                }
                            }
                        }
                    }
                }

                if (nTgt != alli)
                {
                    FatalErrorInFunction << "nTgt:" << nTgt
                        << " alli:" << alli << exit(FatalError);
                }
            }

            // Offsets for slots into results of tgtToSrcMap
            labelList tgtStartOfLocal;
            List<Map<label>> tgtCompactMaps;

            tgtToSrcMap.reset
            (
                new mapDistribute
                (
                    tgtMaps,
                    localRanks,
                    coarseComm,
                    newToOldRanks,
                    tgtStartOfLocal,
                    tgtCompactMaps
                )
            );

            // Assemble srcAddress
            srcAddress.setSize(nSrc);
            if (srcAddress.size())
            {
                label alli = 0;
                forAll(allInterfaces, inti)
                {
                    if (allInterfaces.set(inti))
                    {
                        const auto& intf =
                            refCast<const cyclicACMIGAMGInterface>
                            (
                                allInterfaces[inti]
                            );
                        const auto& AMI = intf.AMI();
                        const auto& srcSlots = AMI.srcAddress();
                        const label localSize =
                            tgtStartOfLocal[inti+1]
                          - tgtStartOfLocal[inti];

                        forAll(srcSlots, srci)
                        {
                            // Append old slots: copy old values and adapt
                            auto& newSlots = srcAddress[alli++];
                            newSlots = srcSlots[srci];
                            // Renumber to new indices
                            mapDistributeBase::renumberMap
                            (
                                newSlots,
                                localSize,
                                tgtStartOfLocal[inti],
                                tgtCompactMaps[inti],
                                AMI.tgtMap().constructHasFlip() //hasFlip
                            );

                            for (const label slot : newSlots)
                            {
                                if
                                (
                                    slot < 0
                                 || slot >= tgtToSrcMap().constructSize()
                                )
                                {
                                    FatalErrorInFunction << " newSlots:"
                                        << newSlots << exit(FatalError);
                                }
                            }
                        }
                    }
                }

                if (nSrc != alli)
                {
                    FatalErrorInFunction << "nSrc:" << nSrc
                        << " alli:" << alli << exit(FatalError);
                }
            }


            // Clean up: if no remote elements sent/received mark as
            // non-distributed. We could do this at the start but this
            // needs to take all the internal transport into account. Easier
            // (but less efficient) to do afterwards now all is compacted.
            {
                const auto& map = srcToTgtMap().subMap();

                bool usesRemote = false;
                forAll(map, proci)
                {
                    if (proci != myAgglom)
                    {
                        const auto& ss = srcToTgtMap().subMap()[proci];
                        const auto& sc = srcToTgtMap().constructMap()[proci];
                        const auto& ts = tgtToSrcMap().subMap()[proci];
                        const auto& tc = tgtToSrcMap().constructMap()[proci];

                        if (ss.size() || sc.size() || ts.size() || tc.size())
                        {
                            usesRemote = true;
                            break;
                        }
                    }
                }

                if (!usesRemote)
                {
                    //Pout<< "** making fully local on new rank "
                    //    << myAgglom << " in comm:" << coarseComm << endl;
                    singlePatchProc = myAgglom;
                    srcToTgtMap.clear();
                    tgtToSrcMap.clear();
                }
            }
        }
        else
        {
            // src/tgt address are straight indices

            srcAddress.setSize(nSrc);
            tgtAddress.setSize(nTgt);

            nSrc = 0;
            nTgt = 0;
            forAll(allInterfaces, inti)
            {
                if (allInterfaces.set(inti))
                {
                    const auto& intf = refCast<const cyclicACMIGAMGInterface>
                    (
                        allInterfaces[inti]
                    );
                    const auto& AMI = intf.AMI();

                    const auto& srcA = AMI.srcAddress();
                    if (srcAddress.size())
                    {
                        label srci = nSrc;
                        forAll(srcA, i)
                        {
                            srcAddress[srci++] = srcA[i]+nTgt;
                        }
                    }

                    const auto& tgtA = AMI.tgtAddress();
                    if (tgtAddress.size())
                    {
                        label tgti = nTgt;
                        forAll(tgtA, i)
                        {
                            tgtAddress[tgti++] = tgtA[i]+nSrc;
                        }
                    }

                    nSrc += srcA.size();
                    nTgt += tgtA.size();
                }
            }
        }

        srcWeights.setSize(nSrc);
        if (hasSrcMagSf)
        {
            srcMagSf.setSize(nSrc);
        }
        if (hasSrcCentroids)
        {
            srcCentroids.setSize(nSrc);
        }
        tgtWeights.setSize(nTgt);
        if (hasTgtMagSf)
        {
            tgtMagSf.setSize(nTgt);
        }


        // Append individual data
        nSrc = 0;
        nTgt = 0;
        forAll(allInterfaces, inti)
        {
            if (allInterfaces.set(inti))
            {
                const auto& intf = refCast<const cyclicACMIGAMGInterface>
                (
                    allInterfaces[inti]
                );
                const auto& AMI = intf.AMI();

                const auto& srcA = AMI.srcAddress();
                {
                    // weights
                    SubList<scalarList>(srcWeights, srcA.size(), nSrc) =
                        AMI.srcWeights();

                    // magSf
                    if (hasSrcMagSf)
                    {
                        SubList<scalar>(srcMagSf, srcA.size(), nSrc) =
                            AMI.srcMagSf();
                    }

                    // centroids
                    if (hasSrcCentroids)
                    {
                        SubList<pointList>(srcCentroids, srcA.size(), nSrc) =
                            AMI.srcCentroids();
                    }
                }

                const auto& tgtA = AMI.tgtAddress();
                {
                    // weights
                    SubList<scalarList>(tgtWeights, tgtA.size(), nTgt) =
                        AMI.tgtWeights();

                    if (hasTgtMagSf)
                    {
                        SubList<scalar>(tgtMagSf, tgtA.size(), nTgt) =
                            AMI.tgtMagSf();
                    }
                }

                nSrc += srcA.size();
                nTgt += tgtA.size();
            }
        }


        // Construct with same arguments as original
        amiPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                AMI.requireMatch(),
                AMI.reverseTarget(),
                AMI.lowWeightCorrection()
            )
        );
        amiPtr_().comm(coarseComm),
        amiPtr_().reset
        (
            std::move(srcToTgtMap),
            std::move(tgtToSrcMap),
            std::move(srcAddress),
            std::move(srcWeights),
            std::move(tgtAddress),
            std::move(tgtWeights),
            singlePatchProc
        );
        amiPtr_().srcMagSf() = std::move(srcMagSf);
        amiPtr_().srcCentroids() = std::move(srcCentroids);
        amiPtr_().tgtMagSf() = std::move(tgtMagSf);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField>
Foam::cyclicACMIGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    const cyclicACMIGAMGInterface& nbr =
        dynamic_cast<const cyclicACMIGAMGInterface&>(neighbPatch());
    const labelUList& nbrFaceCells = nbr.faceCells();

    auto tpnf = tmp<labelField>::New(nbrFaceCells.size());
    labelField& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }

    return tpnf;
}


void Foam::cyclicACMIGAMGInterface::write(Ostream& os) const
{
    GAMGInterface::write(os);

    const bool hasAMI = bool(amiPtr_);

    os  << token::SPACE << neighbPatchID_
        << token::SPACE << owner_
        << token::SPACE << forwardT_
        << token::SPACE << reverseT_
        << token::SPACE << hasAMI;

    if (hasAMI)
    {
        os  << token::SPACE;
        AMI().writeData(os);
    }
}


// ************************************************************************* //
