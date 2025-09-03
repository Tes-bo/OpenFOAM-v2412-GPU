/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ListOps.H"
#include "parLagrangianDistributor.H"
#include "passivePositionParticleCloud.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parLagrangianDistributor::verbose_ = 1;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parLagrangianDistributor::parLagrangianDistributor
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const label nSrcCells,
    const mapDistributePolyMesh& distMap
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap)
{
    const mapDistribute& cellMap = distMap_.cellMap();

    // Get destination processors and cells
    destinationProcID_ = labelList(tgtMesh_.nCells(), UPstream::myProcNo());
    cellMap.reverseDistribute(nSrcCells, destinationProcID_);

    destinationCell_ = identity(tgtMesh_.nCells());
    cellMap.reverseDistribute(nSrcCells, destinationCell_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find all clouds (on all processors) and for each cloud all the objects.
// Result will be synchronised on all processors
void Foam::parLagrangianDistributor::findClouds
(
    const fvMesh& mesh,
    wordList& cloudNames,       // All cloud names on any processor
    boolList& haveClouds,       // Per cloud name, whether my processor has it
    List<wordList>& objectNames // Per cloud name, the field name
)
{
    const IOobject io
    (
        cloud::prefix,
        mesh.time().timeName(),
        mesh.thisDb(),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    // Using the fileHandler:
    // - use fileHandler to synthesise correct processor directory
    // - cannot use readObjects since assumes all processors have same
    //   files (i.e. it only checks processor0)
    const fileNameList localCloudDirs
    (
        fileHandler().readDir
        (
            fileHandler().objectPath
            (
                io,
                word::null  // typeName: not used currently
            ),
            fileName::DIRECTORY
        )
    );


    // Copy out the local cloud names (and fileName -> word)
    cloudNames.resize_nocopy(localCloudDirs.size());
    forAll(localCloudDirs, i)
    {
        cloudNames[i] = localCloudDirs[i];
    }

    // Synchronise cloud names
    Pstream::combineReduce(cloudNames, ListOps::uniqueEqOp<word>());
    Foam::sort(cloudNames);  // Consistent order

    const label nClouds = cloudNames.size();

    // See which of the global cloudNames I have
    haveClouds.resize_nocopy(nClouds);
    haveClouds = false;

    for (const fileName& localCloudName : localCloudDirs)
    {
        const label cloudi = cloudNames.find(localCloudName);
        if (cloudi >= 0)
        {
            haveClouds[cloudi] = true;
        }
    }

    // Collect fields per cloud
    objectNames.resize_nocopy(nClouds);

    for (label cloudi = 0; cloudi < nClouds; ++cloudi)
    {
        objectNames[cloudi].clear();

        if (!haveClouds[cloudi]) continue;

        // Do local scan for valid cloud objects
        const bool oldParRun = UPstream::parRun(false);
        IOobjectList localObjs
        (
            mesh,
            mesh.time().timeName(),
            cloud::prefix/cloudNames[cloudi]
        );
        UPstream::parRun(oldParRun);

        bool isCloud = false;
        if (localObjs.erase("coordinates"))
        {
            isCloud = true;
        }
        if (localObjs.erase("positions"))
        {
            isCloud = true;
        }

        if (isCloud)
        {
            // Has coordinates/positions - so must be a valid cloud
            objectNames[cloudi] = localObjs.sortedNames();
        }
    }

    // Synchronise objectNames (per cloud)
    for (wordList& objNames : objectNames)
    {
        Pstream::combineReduce(objNames, ListOps::uniqueEqOp<word>());
        Foam::sort(objNames);  // Consistent order
    }
}


Foam::autoPtr<Foam::mapDistributeBase>
Foam::parLagrangianDistributor::distributeLagrangianPositions
(
    passivePositionParticleCloud& lpi
) const
{
    //Debug(lpi.size());
    //const label oldLpi = lpi.size();

    labelListList sendMap;

    // Transfer buffers
    PstreamBuffers pBufs;

    {
        // List of lists of particles to be transferred for all of the
        // neighbour processors
        List<IDLList<passivePositionParticle>> particleTransferLists
        (
            UPstream::nProcs()
        );

        // Per particle the destination processor
        labelList destProc(lpi.size());

        label particleI = 0;
        for (passivePositionParticle& ppi : lpi)
        {
            const label destProcI = destinationProcID_[ppi.cell()];
            const label destCellI = destinationCell_[ppi.cell()];

            ppi.cell() = destCellI;
            destProc[particleI++] = destProcI;
            particleTransferLists[destProcI].append(lpi.remove(&ppi));
        }


        // Per processor the indices of the particles to send
        sendMap = invertOneToMany(UPstream::nProcs(), destProc);


        // Stream into send buffers
        forAll(particleTransferLists, procI)
        {
            //Pout<< "To proc " << procI << " sending "
            //    << particleTransferLists[procI] << endl;
            if (particleTransferLists[procI].size())
            {
                UOPstream particleStream(procI, pBufs);
                particleStream << particleTransferLists[procI];
            }
        }
    }


    // Start sending
    pBufs.finishedSends();


    // The cloud name
    const word cloudName = lpi.name();

    {
        // Temporarily rename original cloud so we can construct a new one
        // (to distribute the positions) without getting a duplicate
        // registration warning
        lpi.rename(cloudName + "_old");

        // New empty cloud on tgtMesh
        passivePositionParticleCloud lagrangianPositions
        (
            tgtMesh_,
            Foam::zero{},
            cloudName
        );

        // Retrieve from receive buffers
        for (const int proci : pBufs.allProcs())
        {
            //Pout<< "Receive from processor" << proci << " : "
            //    << pBufs.recvDataCount(proci) << endl;

            if (pBufs.recvDataCount(proci))
            {
                UIPstream particleStream(proci, pBufs);

                // Receive particles and locate them
                IDLList<passivePositionParticle> newParticles
                (
                    particleStream,
                    passivePositionParticle::iNew(tgtMesh_)
                );

                for (passivePositionParticle& newp : newParticles)
                {
                    lagrangianPositions.addParticle(newParticles.remove(&newp));
                }
            }
        }

        const bool writeOnProc = lagrangianPositions.size();

        //if (writeOnProc)
        {
            // Write coordinates file
            IOPosition<passivePositionParticleCloud>
            (
                lagrangianPositions
            ).write(writeOnProc);

            // Optionally write positions file in v1706 format and earlier
            if (particle::writeLagrangianPositions)
            {
                IOPosition<passivePositionParticleCloud>
                (
                     lagrangianPositions,
                     cloud::geometryType::POSITIONS
                 ).write(writeOnProc);
            }
        }
        //else if (!writeOnProc && oldLpi)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    // Remove any existing coordinates
        //    Foam::rm
        //    (
        //        IOPosition<passivePositionParticleCloud>
        //        (
        //            lagrangianPositions
        //        ).objectPath()
        //    );
        //
        //    // Remove any existing positions
        //    Foam::rm
        //    (
        //        IOPosition<passivePositionParticleCloud>
        //        (
        //            lagrangianPositions,
        //            cloud::geometryType::POSITIONS
        //        ).objectPath()
        //    );
        //}
    }

    // Restore cloud name
    lpi.rename(cloudName);


    // The constructMap is in linear (processor) order
    return autoPtr<mapDistributeBase>::New
    (
        mapDistributeBase::layoutTypes::linear,
        std::move(sendMap)
    );
}


Foam::autoPtr<Foam::mapDistributeBase>
Foam::parLagrangianDistributor::distributeLagrangianPositions
(
    const word& cloudName
) const
{
    // Mixed exists/missing on various ranks?
    // Avoid masterRead+broadcast (can cause blocking)

    auto& handler = Foam::fileHandler();
    const bool oldDistributed =
        handler.distributed
        (
            !fileOperation::cacheLevel() || handler.distributed()
        );


    // Load cloud
    passivePositionParticleCloud lpi(srcMesh_, cloudName, false);

    // Restore distributed flag
    handler.distributed(oldDistributed);

    // Distribute particles to other ranks
    return distributeLagrangianPositions(lpi);
}


// ************************************************************************* //
