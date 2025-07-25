/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "loadOrCreateMesh.H"
#include "faMesh.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "decomposedBlockData.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::checkFileExistence(const fileName& fName)
{
    // Trimmed-down version of lookupAndCacheProcessorsPath
    // with Foam::exists() check. No caching.

    // Check for two conditions:
    // - file has to exist
    // - if collated the entry has to exist inside the file

    // Note: bypass fileOperation::filePath(IOobject&) since has problems
    //       with going to a different number of processors
    //       (in collated format). Use file-based searching instead

    const auto& handler = Foam::fileHandler();
    typedef fileOperation::procRangeType procRangeType;

    fileName path, pDir, local;
    procRangeType group;
    label numProcs;
    const label proci =
        fileOperation::splitProcessorPath
        (fName, path, pDir, local, group, numProcs);

    bool found = false;

    if (proci != -1)
    {
        // Read all directories to see any beginning with processor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const fileNameList dirEntries
        (
            handler.readDir(path, fileName::Type::DIRECTORY)
        );

        // Extract info from processorN or processorsNN
        // - highest processor number
        // - directory+offset containing data for proci

        // label nProcs = 0;
        for (const fileName& dirN : dirEntries)
        {
            // Analyse directory name
            fileName rp, rd, rl;
            label rNum;
            const label readProci =
                fileOperation::splitProcessorPath
                (dirN, rp, rd, rl, group, rNum);

            if (proci == readProci)
            {
                // Found "processorN"
                if (Foam::exists(path/dirN/local))
                {
                    found = true;
                    break;
                }
            }
            else if (rNum != -1)
            {
                // "processorsNN" or "processorsNN_start-end"
                if (group.empty())
                {
                    // "processorsNN"
                    if (proci < rNum && Foam::exists(path/dirN/local))
                    {
                        found = true;
                        break;
                    }
                }
                else if (group.contains(proci))
                {
                    // "processorsNN_start-end"
                    // - save the local proc offset

                    if (Foam::exists(path/dirN/local))
                    {
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    if (!found)
    {
        found = Foam::exists(fName);
    }

    return found;
}


Foam::boolList Foam::haveMeshFile
(
    const Time& runTime,
    const fileName& meshPath,
    const word& meshFile,
    const bool verbose
)
{
    #if 0

    // Simple directory scanning - too fragile
    bool found = checkFileExistence(runTime.path()/meshPath/meshFile);

    #else

    // Trimmed-down version of lookupAndCacheProcessorsPath
    // with Foam::exists() check. No caching.

    // Check for two conditions:
    // - file has to exist
    // - if collated the entry has to exist inside the file

    // Note: bypass fileOperation::filePath(IOobject&) since has problems
    //       with going to a different number of processors
    //       (in collated format). Use file-based searching instead

    const auto& handler = Foam::fileHandler();
    typedef fileOperation::procRangeType procRangeType;

    const fileName fName
    (
        handler.filePath(runTime.path()/meshPath/meshFile)
    );
    bool found = handler.isFile(fName);
    if (returnReduceAnd(found)) // worldComm
    {
        // Bit tricky: avoid having all slaves open file since this involves
        // reading it on master and broadcasting it. This fails if file > 2G.
        // So instead only read on master

        bool isCollated = false;

        // Note: can test only world-master. Since even host-collated will have
        // same file format type for all processors
        if (UPstream::master(UPstream::worldComm))
        {
            const bool oldParRun = UPstream::parRun(false);

            IFstream is(fName);
            if (is.good())
            {
                IOobject io(meshFile, meshPath, runTime);
                io.readHeader(is);

                isCollated = decomposedBlockData::isCollatedType(io);
            }
            UPstream::parRun(oldParRun);
        }
        Pstream::broadcast(isCollated); //UPstream::worldComm


        // Collect block-number in individual filenames (might differ
        // on different processors)
        if (isCollated)
        {
            const label nProcs = UPstream::nProcs(fileHandler().comm());
            const label myProcNo = UPstream::myProcNo(fileHandler().comm());

            // Collect file names on master of local communicator
            const fileNameList fNames
            (
                Pstream::listGatherValues
                (
                    fName,
                    fileHandler().comm(),
                    UPstream::msgType()
                )
            );

            // Collect local block number
            label myBlockNumber = -1;
            {
                fileName path, pDir, local;
                procRangeType group;
                label numProcs;
                label proci = fileOperation::splitProcessorPath
                (
                    fName,
                    path,
                    pDir,
                    local,
                    group,
                    numProcs
                );

                if (proci == -1 && group.empty())
                {
                    // 'processorsXXX' format so contains all ranks
                    // according to worldComm
                    myBlockNumber = UPstream::myProcNo(UPstream::worldComm);
                }
                else
                {
                    // 'processorsXXX_n-m' format so check for the
                    // relative rank
                    myBlockNumber = myProcNo;
                }
            }
            const labelList myBlockNumbers
            (
                Pstream::listGatherValues
                (
                    myBlockNumber,
                    fileHandler().comm(),
                    UPstream::msgType()
                )
            );



            // Determine for all whether the filename exists in the collated
            // file.
            boolList allFound(nProcs, false);

            if (UPstream::master(fileHandler().comm()))
            {
                // Store nBlocks and index of file that was used for nBlocks
                label nBlocks = -1;
                label blockRanki = -1;
                forAll(fNames, ranki)
                {
                    if
                    (
                        blockRanki == -1
                     || (fNames[ranki] != fNames[blockRanki])
                    )
                    {
                        blockRanki = ranki;
                        IFstream is(fNames[ranki]);
                        nBlocks = decomposedBlockData::getNumBlocks(is);
                    }

                    allFound[ranki] = (myBlockNumbers[ranki] < nBlocks);
                }
            }

            found = Pstream::listScatterValues
            (
                allFound,
                fileHandler().comm(),
                UPstream::msgType()
            );
        }
    }
    #endif

    // Globally consistent information about who has a mesh
    boolList haveFileOnProc
    (
        UPstream::allGatherValues<bool>(found, UPstream::worldComm)
    );

    if (verbose)
    {
        Info<< "Per processor availability of \""
            << meshFile << "\" file in " << meshPath << nl
            << "    " << flatOutput(haveFileOnProc) << nl << endl;
    }

    return haveFileOnProc;
}


void Foam::removeProcAddressing(const faMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "edge", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::removeProcAddressing(const polyMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "cell", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::masterMeshInstance
(
    const IOobject& io,
    fileName& facesInstance,
    fileName& pointsInstance
)
{
    const fileName meshSubDir
    (
        polyMesh::meshDir(io.name())
    );

    if (UPstream::master())
    {
        const bool oldParRun = UPstream::parRun(false);
        const label oldNumProcs = fileHandler().nProcs();
        const int oldCache = fileOperation::cacheLevel(0);

        facesInstance = io.time().findInstance
        (
            meshSubDir,
            "faces",
            IOobjectOption::MUST_READ
        );
        pointsInstance = io.time().findInstance
        (
            meshSubDir,
            "points",
            IOobjectOption::MUST_READ
        );

        fileOperation::cacheLevel(oldCache);
        if (oldParRun)
        {
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
        }
        UPstream::parRun(oldParRun);
    }

    // Broadcast information to all
    Pstream::broadcasts
    (
        UPstream::worldComm,
        facesInstance,
        pointsInstance
    );
}


// ************************************************************************* //
