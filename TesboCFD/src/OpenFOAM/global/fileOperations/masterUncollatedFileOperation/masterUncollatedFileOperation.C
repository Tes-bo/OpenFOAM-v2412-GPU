/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "masterUncollatedFileOperation.H"
#include "fileOperationInitialise.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "instant.H"
#include "IFstream.H"
#include "SpanStream.H"
#include "masterOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "dummyISstream.H"
#include "SubList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(masterUncollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        masterUncollatedFileOperation,
        word
    );
    addToRunTimeSelectionTable
    (
        fileOperation,
        masterUncollatedFileOperation,
        comm
    );

    float masterUncollatedFileOperation::maxMasterFileBufferSize
    (
        Foam::debug::floatOptimisationSwitch("maxMasterFileBufferSize", 1e9)
    );
    registerOptSwitch
    (
        "maxMasterFileBufferSize",
        float,
        masterUncollatedFileOperation::maxMasterFileBufferSize
    );

    // Threaded MPI: not required
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        fileOperationInitialise_unthreaded,
        word,
        masterUncollated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::filePathInfo
(
    const bool checkGlobal,
    const bool isFile,
    const IOobject& io,
    const dirIndexList& pDirs,
    const bool search,
    pathType& searchType,
    word& procsDir,
    word& newInstancePath
) const
{
    procsDir.clear();
    newInstancePath.clear();

    if (io.instance().isAbsolute())
    {
        fileName objPath = io.instance()/io.name();

        if (isFileOrDir(isFile, objPath))
        {
            searchType = fileOperation::ABSOLUTE;
            return objPath;
        }
        else
        {
            searchType = fileOperation::NOTFOUND;
            return fileName();
        }
    }
    else
    {
        // 1. Check the writing fileName
        fileName writePath(objectPath(io, io.headerClassName()));

        if (isFileOrDir(isFile, writePath))
        {
            searchType = fileOperation::WRITEOBJECT;
            return writePath;
        }

        // 2. Check processors/
        if (io.time().processorCase())
        {
            for (const dirIndex& dirIdx : pDirs)
            {
                const fileName& pDir = dirIdx.first();
                fileName objPath =
                    processorsPath(io, io.instance(), pDir)
                   /io.name();
                if (objPath != writePath && isFileOrDir(isFile, objPath))
                {
                    searchType = dirIdx.second().first();
                    procsDir = pDir;
                    return objPath;
                }
            }
        }
        {
            // 3. Check local
            fileName localPath = io.objectPath();

            if
            (
                localPath != writePath
            &&  isFileOrDir(isFile, localPath)
            )
            {
                searchType = fileOperation::OBJECT;
                return localPath;
            }
        }



        // Any global checks
        if
        (
            checkGlobal
         && io.time().processorCase()
         && (
                io.instance() == io.time().system()
             || io.instance() == io.time().constant()
            )
        )
        {
            fileName parentPath =
            (
                io.rootPath()/io.globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name()
            );

            if (isFileOrDir(isFile, parentPath))
            {
                searchType = fileOperation::PARENTOBJECT;
                return parentPath;
            }
        }

        // Check for approximately same time. E.g. if time = 1e-2 and
        // directory is 0.01 (due to different time formats)
        const auto pathFnd = times_.cfind(io.time().path());

        if (search && pathFnd.good())
        {
            newInstancePath =
                Time::findInstancePath
                (
                    *pathFnd(),
                    instant(io.instance())
                );

            if (newInstancePath.size() && newInstancePath != io.instance())
            {
                // 1. Try processors equivalent
                for (const dirIndex& dirIdx : pDirs)
                {
                    const fileName& pDir = dirIdx.first();

                    fileName fName
                    (
                        processorsPath(io, newInstancePath, pDir)
                       /io.name()
                    );
                    if (isFileOrDir(isFile, fName))
                    {
                        switch (dirIdx.second().first())
                        {
                            case fileOperation::PROCUNCOLLATED:
                            {
                                searchType =
                                    fileOperation::PROCUNCOLLATEDINSTANCE;
                            }
                            break;
                            case fileOperation::PROCBASEOBJECT:
                            {
                                searchType = fileOperation::PROCBASEINSTANCE;
                            }
                            break;
                            case fileOperation::PROCOBJECT:
                            {
                                searchType = fileOperation::PROCINSTANCE;
                            }
                            break;
                            default:
                            break;
                        }
                        procsDir = pDir;
                        return fName;
                    }
                }


                // 2. Check local
                fileName fName
                (
                   io.rootPath()/io.caseName()
                  /newInstancePath/io.db().dbDir()/io.local()/io.name()
                );
                if (isFileOrDir(isFile, fName))
                {
                    searchType = fileOperation::FINDINSTANCE;
                    return fName;
                }
            }
        }
    }

    // Nothing found
    searchType = fileOperation::NOTFOUND;
    return fileName();
}


Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::localObjectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& procDir,
    const word& instancePath
) const
{
    // Replacement for IOobject::objectPath()

    switch (searchType)
    {
        case fileOperation::ABSOLUTE:
        {
            return io.instance()/io.name();
        }
        break;

        case fileOperation::OBJECT:
        {
            return io.path()/io.name();
        }
        break;

        case fileOperation::WRITEOBJECT:
        {
            return objectPath(io, io.headerClassName());
        }
        break;

        case fileOperation::PROCUNCOLLATED:
        {
            // Uncollated type, e.g. processor1
            const word procName
            (
                "processor" + Foam::name(Pstream::myProcNo(UPstream::worldComm))
            );
            return
                processorsPath
                (
                    io,
                    io.instance(),
                    (
                        Pstream::parRun()
                      ? procName
                      : procDir
                    )
                )
               /io.name();
        }
        break;

        case fileOperation::PROCBASEOBJECT:
        {
            // Collated, e.g. processors4
            return
                processorsPath(io, io.instance(), procDir)
               /io.name();
        }
        break;

        case fileOperation::PROCOBJECT:
        {
            // Processors directory locally provided by the fileHandler itself
            return
                processorsPath(io, io.instance(), processorsDir(io))
               /io.name();
        }
        break;

        case fileOperation::PARENTOBJECT:
        {
            return
                io.rootPath()/io.globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::FINDINSTANCE:
        {
            return
                io.rootPath()/io.caseName()
               /instancePath/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::PROCUNCOLLATEDINSTANCE:
        {
            // Uncollated type, e.g. processor1
            const word procName
            (
                "processor"
              + Foam::name(Pstream::myProcNo(UPstream::worldComm))
            );
            return
                processorsPath
                (
                    io,
                    instancePath,
                    (
                        Pstream::parRun()
                      ? procName
                      : procDir
                    )
                )
               /io.name();
        }
        break;

        case fileOperation::PROCBASEINSTANCE:
        {
            // Collated, e.g. processors4
            return
                processorsPath(io, instancePath, procDir)
               /io.name();
        }
        break;

        case fileOperation::PROCINSTANCE:
        {
            // Processors directory locally provided by the fileHandler itself
            return
                processorsPath(io, instancePath, processorsDir(io))
               /io.name();
        }
        break;

        case fileOperation::NOTFOUND:
        {
            return fileName();
        }
        break;

        default:
        {
            NotImplemented;
            return fileName();
        }
    }
}


void Foam::fileOperations::masterUncollatedFileOperation::readAndSend
(
    const fileName& filePath,
    const labelUList& recvProcs,
    PstreamBuffers& pBufs
)
{
    if (recvProcs.empty()) return;

    IFstream ifs(filePath, IOstreamOption::BINARY);

    if (!ifs.good())
    {
        FatalIOErrorInFunction(filePath)
            << "Cannot open file " << filePath
            << exit(FatalIOError);
    }

    // Read file contents (compressed or uncompressed) into a character buffer
    DynamicList<char> buf(IFstream::readContents(ifs));

    if (debug)
    {
        Info<< "masterUncollatedFileOperation::readAndSend :"
            << filePath
            << " (compressed:" << bool(ifs.compression())
            << ") : " << " bytes" << endl;
    }

    for (const label proci : recvProcs)
    {
        UOPstream os(proci, pBufs);
        os.write(buf.cdata_bytes(), buf.size_bytes());
    }
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::read
(
    IOobject& io,
    const label comm,
    const bool uniform,             // on comms master only
    const fileNameList& filePaths,  // on comms master and sub-ranks
    const boolUList& readOnProcs    // on comms master and sub-ranks
)
{
    autoPtr<ISstream> isPtr;

    PstreamBuffers pBufs(comm);

    if (UPstream::master(comm))
    {
        if (uniform)
        {
            if (readOnProcs[0])
            {
                if (filePaths[0].empty())
                {
                    FatalIOErrorInFunction(filePaths[0])
                        << "Cannot find file " << io.objectPath()
                        << " fileHandler : comm:" << comm
                        << " ioRanks:" << UPstream::procID(comm)
                        << exit(FatalIOError);
                }

                DynamicList<label> recvProcs(UPstream::nProcs(comm));
                for (const int proci : UPstream::allProcs(comm))
                {
                    if (readOnProcs[proci])
                    {
                        recvProcs.push_back(proci);
                    }
                }

                // Read on master and send to all processors
                // (including master for simplicity)
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream :"
                        << " For uniform file " << filePaths[0]
                        << " sending to " << recvProcs
                        << " in comm:" << comm << endl;
                }
                readAndSend(filePaths[0], recvProcs, pBufs);
            }
        }
        else
        {
            if (readOnProcs[0])
            {
                if (filePaths[0].empty())
                {
                    FatalIOErrorInFunction(filePaths[0])
                        << "Cannot find file " << io.objectPath()
                        << " fileHandler : comm:" << comm
                        << " ioRanks:" << UPstream::procID(comm)
                        << exit(FatalIOError);
                }

                // Open master
                isPtr.reset(new IFstream(filePaths[0]));

                // Read header
                if (!io.readHeader(*isPtr))
                {
                    FatalIOErrorInFunction(*isPtr)
                        << "problem while reading header for object "
                        << io.name()
                        << " fileHandler : comm:" << comm
                        << " ioRanks:" << UPstream::procID(comm)
                        << exit(FatalIOError);
                }
            }

            // Read sub-rank files
            for (const int proci : UPstream::subProcs(comm))
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream :"
                        << " For processor " << proci
                        << " opening " << filePaths[proci] << endl;
                }

                const fileName& fPath = filePaths[proci];

                if (readOnProcs[proci] && !fPath.empty())
                {
                    // Note: handle compression ourselves since size cannot
                    // be determined without actually uncompressing
                    readAndSend(fPath, labelList(one{}, proci), pBufs);
                }
            }
        }
    }

    pBufs.finishedScatters();

    // isPtr will be valid on master and will be the unbuffered
    // IFstream. Else the information is in the PstreamBuffers (and
    // the special case of a uniform file)

    if (!isPtr)
    {
        if (readOnProcs[UPstream::myProcNo(comm)])
        {
            // This processor needs to return something
            List<char> buf(pBufs.recvDataCount(UPstream::masterNo()));

            if (!buf.empty())
            {
                UIPstream is(UPstream::masterNo(), pBufs);
                is.read(buf.data(), buf.size());
            }

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::readStream :"
                    << " Done reading " << buf.size() << " bytes" << endl;
            }

            // A local character buffer copy of the Pstream contents.
            // Construct with same parameters (ASCII, current version)
            // as the IFstream so that it has the same characteristics.

            isPtr.reset(new ICharStream(std::move(buf)));

            // With the proper file name
            isPtr->name() = filePaths[UPstream::myProcNo(comm)];

            if (!io.readHeader(*isPtr))
            {
                FatalIOErrorInFunction(*isPtr)
                    << "problem while reading header for object "
                    << io.name()
                    << " fileHandler : comm:" << comm
                    << " ioRanks:" << UPstream::procID(comm)
                    << exit(FatalIOError);
            }
        }
        else
        {
            isPtr.reset(new dummyISstream());
        }
    }

    return isPtr;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Construction helper: self/world/local communicator and IO ranks
static Tuple2<label, labelList> getCommPattern()
{
    // Default is COMM_WORLD (single master)
    Tuple2<label, labelList> commAndIORanks
    (
        UPstream::worldComm,
        fileOperation::getGlobalIORanks()
    );

    if (UPstream::parRun() && commAndIORanks.second().size() > 1)
    {
        // Multiple masters: ranks for my IO range
        commAndIORanks.first() = UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            fileOperation::subRanks(commAndIORanks.second())
        );
    }

    return commAndIORanks;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::fileOperations::masterUncollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        DetailInfo
            << "I/O    : " << typeName
            << " (maxMasterFileBufferSize " << maxMasterFileBufferSize << ')'
            << endl;
    }

    if (IOobject::fileModificationChecking == IOobject::timeStampMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
        IOobject::fileModificationChecking = IOobject::timeStamp;
    }
    else if (IOobject::fileModificationChecking == IOobject::inotifyMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify"
                << endl;
        }
        IOobject::fileModificationChecking = IOobject::inotify;
    }
}


Foam::fileOperations::masterUncollatedFileOperation::
masterUncollatedFileOperation
(
    bool verbose
)
:
    fileOperation
    (
        getCommPattern()
    ),
    managedComm_(getManagedComm(comm_))  // Possibly locally allocated
{
    init(verbose);

    if (comm_ == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_ << exit(FatalError);
    }
    if (UPstream::nProcs(comm_) == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_
            << " nProcs:" << UPstream::nProcs(comm_)
            << exit(FatalError);
    }
    if (UPstream::myProcNo(comm_) == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_
            << " myProcNo:" << UPstream::myProcNo(comm_)
            << exit(FatalError);
    }
}


Foam::fileOperations::masterUncollatedFileOperation::
masterUncollatedFileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
:
    fileOperation(commAndIORanks, distributedRoots),
    managedComm_(-1)  // Externally managed
{
    init(verbose);

    if (comm_ == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_ << exit(FatalError);
    }
    if (UPstream::nProcs(comm_) == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_
            << " nProcs:" << UPstream::nProcs(comm_)
            << exit(FatalError);
    }
    if (UPstream::myProcNo(comm_) == -1)
    {
        FatalErrorInFunction<< "Problem comm_:" << comm_
            << " myProcNo:" << UPstream::myProcNo(comm_)
            << exit(FatalError);
    }
}


void Foam::fileOperations::masterUncollatedFileOperation::storeComm() const
{
    // From externally -> locally managed
    managedComm_ = getManagedComm(comm_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterUncollatedFileOperation::
~masterUncollatedFileOperation()
{
    UPstream::freeCommunicator(managedComm_);
}


// * * * * * * * * * * * * * Filesystem Operations * * * * * * * * * * * * * //

bool Foam::fileOperations::masterUncollatedFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterOp<bool>
    (
        dir,
        mkDirOp(mode),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterOp<bool>
    (
        fName,
        chModOp(mode),
        UPstream::msgType(),
        comm_
    );
}


mode_t Foam::fileOperations::masterUncollatedFileOperation::mode
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<mode_t>
    (
        fName,
        modeOp(followLink),
        UPstream::msgType(),
        comm_
    );
}


Foam::fileName::Type Foam::fileOperations::masterUncollatedFileOperation::type
(
    const fileName& fName,
    const bool followLink
) const
{
    return fileName::Type
    (
        masterOp<int>
        (
            fName,
            typeOp(followLink),
            UPstream::msgType(),
            comm_
        )
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return masterOp<bool>
    (
        fName,
        existsOp(checkGzip, followLink),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<bool>
    (
        fName,
        isDirOp(followLink),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return masterOp<bool>
    (
        fName,
        isFileOp(checkGzip, followLink),
        UPstream::msgType(),
        comm_
    );
}


off_t Foam::fileOperations::masterUncollatedFileOperation::fileSize
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<off_t>
    (
        fName,
        fileSizeOp(followLink),
        UPstream::msgType(),
        comm_
    );
}


time_t Foam::fileOperations::masterUncollatedFileOperation::lastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<time_t>
    (
        fName,
        lastModifiedOp(followLink),
        UPstream::msgType(),
        UPstream::worldComm
    );
}


double Foam::fileOperations::masterUncollatedFileOperation::highResLastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<double>
    (
        fName,
        highResLastModifiedOp(followLink),
        UPstream::msgType(),
        UPstream::worldComm
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterOp<bool>
    (
        fName,
        mvBakOp(ext),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::rm
(
    const fileName& fName
) const
{
    return masterOp<bool>
    (
        fName,
        rmOp(),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::rmDir
(
    const fileName& dir,
    const bool silent,
    const bool emptyOnly
) const
{
    return masterOp<bool>
    (
        dir,
        rmDirOp(silent, emptyOnly),
        UPstream::msgType(),
        comm_
    );
}


Foam::fileNameList Foam::fileOperations::masterUncollatedFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
) const
{
    return masterOp<fileNameList>
    (
        dir,
        readDirOp(type, filtergz, followLink),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool>
    (
        src,
        dst,
        cpOp(followLink),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool>
    (
        src,
        dst,
        lnOp(),
        UPstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool>
    (
        src,
        dst,
        mvOp(followLink),
        UPstream::msgType(),
        comm_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName,
    const bool search
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal
            << " parRun:" << Pstream::parRun()
            << " localmaster:" << Pstream::master(comm_) << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    const refPtr<dirIndexList> pDirs(lookupProcessorsPath(io.objectPath()));

    // Trigger caching of times
    if (cacheLevel() > 0)
    {
        (void)findTimes(io.time().path(), io.time().constant());
    }

    // Determine master filePath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    if (Pstream::master(comm_))
    {
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        // All masters search locally. Note that global objects might
        // fail (except on master). This gets handled later on (in PARENTOBJECT)
        objPath =
            filePathInfo
            (
                checkGlobal,
                true,
                io,
                pDirs,
                search,
                searchType,
                procsDir,
                newInstancePath
            );

        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::filePath :"
                << " master objPath:" << objPath
                << " searchType:" << fileOperation::pathTypeNames_[searchType]
                << " procsDir:" << procsDir << " instance:" << newInstancePath
                << endl;
        }
    }

    // Broadcast information about where the master found the object
    // Note: use the worldComm to make sure all processors decide
    //       the same type. Only procsDir is allowed to differ; searchType
    //       and instance have to be same
    if (UPstream::parRun())
    {
        int masterType(searchType);
        Pstream::broadcasts(UPstream::worldComm, masterType, newInstancePath);
        searchType = pathType(masterType);
    }

    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
        // Distribute master path. This makes sure it is seen as uniform
        // and only gets read from the master.
        Pstream::broadcasts(UPstream::worldComm, objPath, procsDir);
    }
    else
    {
        Pstream::broadcast(procsDir, comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)

                objPath = masterOp<fileName>
                (
                    io.objectPath(),
                    fileOrNullOp(true), // isFile=true
                    UPstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::filePath :"
            << " Returning from file searching using type "
            << fileOperation::pathTypeNames_[searchType] << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io,
    const bool search
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::dirPath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal
            << " parRun:" << Pstream::parRun()
            << " localmaster:" << Pstream::master(comm_) << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    const refPtr<dirIndexList> pDirs(lookupProcessorsPath(io.objectPath()));

    // Trigger caching of times
    if (cacheLevel() > 0)
    {
        (void)findTimes(io.time().path(), io.time().constant());
    }

    // Determine master dirPath and broadcast

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    // Local IO node searches for file
    if (Pstream::master(comm_))
    {
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        objPath = filePathInfo
        (
            checkGlobal,
            false,
            io,
            pDirs,
            search,
            searchType,
            procsDir,
            newInstancePath
        );

        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::dirPath :"
                << " master objPath:" << objPath
                << " searchType:" << fileOperation::pathTypeNames_[searchType]
                << " procsDir:" << procsDir << " instance:" << newInstancePath
                << endl;
        }
    }


    // Broadcast information about where the master found the object
    // Note: use the worldComm to make sure all processors decide
    //       the same type. Only procsDir is allowed to differ; searchType
    //       and instance have to be same
    if (UPstream::parRun())
    {
        int masterType(searchType);
        Pstream::broadcasts(UPstream::worldComm, masterType, newInstancePath);
        searchType = pathType(masterType);
    }


    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
        // Distribute master path. This makes sure it is seen as uniform
        // and only gets read from the master.
        Pstream::broadcasts(UPstream::worldComm, objPath, procsDir);
    }
    else
    {
        // Broadcast local processors dir amongst all local nodes
        Pstream::broadcast(procsDir, comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)

                objPath = masterOp<fileName>
                (
                    io.objectPath(),
                    fileOrNullOp(false), // isFile=false
                    UPstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::dirPath :"
            << " Returning from directory searching using type "
            << fileOperation::pathTypeNames_[searchType] << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::masterUncollatedFileOperation::exists
(
    const dirIndexList& pDirs,
    IOobject& io
) const
{
    // Cut-down version of filePathInfo that does not look for
    // different instance or parent directory

    const bool isFile = !io.name().empty();

    // Generate output filename for object
    const fileName writePath(objectPath(io, word::null));

    // 1. Test writing name for either directory or a (valid) file
    if (isFileOrDir(isFile, writePath))
    {
        return true;
    }

    // 2. Check processors/
    if (io.time().processorCase())
    {
        for (const dirIndex& dirIdx : pDirs)
        {
            const fileName& pDir = dirIdx.first();
            fileName procPath =
                processorsPath(io, io.instance(), pDir)
               /io.name();
            if (procPath != writePath && isFileOrDir(isFile, procPath))
            {
                return true;
            }
        }
    }

    // 3. Check local
    fileName localPath = io.objectPath();

    if (localPath != writePath && isFileOrDir(isFile, localPath))
    {
        return true;
    }

    return false;
}


Foam::IOobject
Foam::fileOperations::masterUncollatedFileOperation::findInstance
(
    const IOobject& startIO,
    const scalar startValue,
    const word& stopInstance,
    const bool constant_fallback
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::findInstance :"
            << " Starting searching for name:" << startIO.name()
            << " local:" << startIO.local()
            << " from instance:" << startIO.instance()
            << endl;
    }

    const Time& time = startIO.time();
    IOobject io(startIO);

    // Note: - if name is empty, just check the directory itself
    //       - check both for isFile and headerOk since the latter does a
    //         filePath so searches for the file.
    //       - check for an object with local file scope (so no looking up in
    //         parent directory in case of parallel)


    const refPtr<dirIndexList> pDirs(lookupProcessorsPath(io.objectPath()));

    word foundInstance;

    // if (Pstream::master(comm_))
    if (Pstream::master(UPstream::worldComm))
    {
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        if (exists(pDirs, io))
        {
            foundInstance = io.instance();
        }
        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);
    }


    // Do parallel early exit to avoid calling time.times()
    Pstream::broadcast(foundInstance, UPstream::worldComm);

    if (!foundInstance.empty())
    {
        io.instance() = foundInstance;
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findInstance :"
                << " for name:" << io.name() << " local:" << io.local()
                << " found starting instance:" << io.instance() << endl;
        }
        return io;
    }


    // Handling failures afterwards
    const bool exitIfMissing = startIO.isReadRequired();

    enum failureCodes { FAILED_STOPINST = 1, FAILED_CONSTINST = 2 };
    int failed(0);

    instantList ts = time.times();

    // if (Pstream::master(comm_))
    if (Pstream::master(UPstream::worldComm))
    {
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        label instIndex = ts.size()-1;

        // Backward search for first time that is <= startValue
        for (; instIndex >= 0; --instIndex)
        {
            if (ts[instIndex].value() <= startValue)
            {
                break;
            }
        }

        // Continue (forward) searching from here
        for (; instIndex >= 0; --instIndex)
        {
            // Shortcut: if actual directory is the timeName we've
            // already tested it
            if (ts[instIndex].name() == time.timeName())
            {
                continue;
            }

            io.instance() = ts[instIndex].name();
            if (exists(pDirs, io))
            {
                foundInstance = io.instance();
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::findInstance :"
                        << " for name:" << io.name() << " local:" << io.local()
                        << " found at:" << io.instance()
                        << endl;
                }
                break;
            }

            // Check if hit minimum instance
            if (io.instance() == stopInstance)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::findInstance :"
                        << " name:" << io.name()
                        << " local:" << io.local()
                        << " at stop-instance:" << io.instance() << endl;
                }

                if (exitIfMissing)
                {
                    failed = failureCodes::FAILED_STOPINST;
                }
                else
                {
                    // At the stopInstance
                    foundInstance = io.instance();
                }
                break;
            }
        }


        // times() usually already includes the constant() so would
        // have been checked above. However, re-test under these conditions:
        // - times() is empty. Sometimes this can happen (e.g. decomposePar
        //   with collated)
        // - times()[0] is not constant
        // - Times is empty.
        //   Sometimes this can happen (eg, decomposePar with collated)
        // - Times[0] is not constant
        // - The startValue is negative (eg, kivaTest).
        //   This plays havoc with the reverse search, causing it to miss
        //   'constant'

        if
        (
            !failed && foundInstance.empty()
         && (ts.empty() || ts[0].name() != time.constant() || startValue < 0)
        )
        {
            // Note. This needs to be a hard-coded "constant" (not constant
            // function of Time), because the latter points to
            // the case constant directory in parallel cases.
            // However, parRun is disabled so they are actually the same.

            io.instance() = time.constant();

            if (exists(pDirs, io))
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::findInstance :"
                        << " name:" << io.name()
                        << " local:" << io.local()
                        << " at:" << io.instance() << endl;
                }
                foundInstance = io.instance();
            }
        }

        if (!failed && foundInstance.empty())
        {
            if (exitIfMissing)
            {
                failed = failureCodes::FAILED_CONSTINST;
            }
            else if (constant_fallback)
            {
                foundInstance = time.constant();
            }
            else
            {
                foundInstance.clear();
            }
        }

        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);  // Restore parallel state
    }

    Pstream::broadcast(foundInstance, UPstream::worldComm);

    io.instance() = foundInstance;


    // Handle failures
    // ~~~~~~~~~~~~~~~
    if (failed)
    {
        FatalErrorInFunction << "Cannot find";

        if (!io.name().empty())
        {
            FatalError
                << " file \"" << io.name() << "\" in";
        }

        FatalError
            << " directory "
            << io.local() << " in times "
            << startIO.instance() << " down to ";

        if (failed == failureCodes::FAILED_STOPINST)
        {
            FatalError << stopInstance;
        }
        else  // FAILED_CONSTINST
        {
            FatalError << "constant";
        }
        FatalError << exit(FatalError);
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::findInstance :"
            << " name:" << io.name() << " local:" << io.local()
            << " returning instance:" << io.instance() << endl;
    }
    return io;
}


Foam::fileNameList
Foam::fileOperations::masterUncollatedFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readObjects :"
            << " db:" << db.objectPath()
            << " local:" << local << " instance:" << instance << endl;
    }

    fileNameList objectNames;
    newInstance.clear();

    // Note: readObjects uses WORLD to make sure order of objects is the
    //       same everywhere

    if (Pstream::master(UPstream::worldComm))
    {
        // Avoid fileOperation::readObjects from triggering parallel ops
        // (through call to filePath which triggers parallel )
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        //- Use non-time searching version
        objectNames = fileOperation::readObjects
        (
            db,
            instance,
            local,
            newInstance
        );

        if (newInstance.empty())
        {
            // Find similar time

            // Copy of Time::findInstancePath. We want to avoid the
            // parallel call to findTimes.
            const instantList timeDirs
            (
                fileOperation::findTimes
                (
                    db.time().path(),
                    db.time().constant()
                )
            );

            fileName foundInst
            (
                Time::findInstancePath(timeDirs, instant(instance))
            );

            if (!foundInst.empty())
            {
                objectNames = fileOperation::readObjects
                (
                    db,
                    foundInst,     // newly found time
                    local,
                    newInstance
                );
            }
        }

        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);  // Restore parallel state
    }

    Pstream::broadcasts(UPstream::worldComm, newInstance, objectNames);

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readObjects :"
            << " newInstance:" << newInstance
            << " objectNames:" << objectNames << endl;
    }

    return objectNames;
}


bool Foam::fileOperations::masterUncollatedFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    bool ok = false;

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readHeader :" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << fName << endl;
    }

    // We assume if filePath is the same
    //  - headerClassName
    //  - note
    // are also the same, independent of where the file came from.

    // Get filePaths on world master
    fileNameList filePaths(Pstream::nProcs(UPstream::worldComm));
    filePaths[UPstream::myProcNo(UPstream::worldComm)] = fName;
    Pstream::gatherList(filePaths, UPstream::msgType(), UPstream::worldComm);

    bool uniform
    (
        UPstream::master(UPstream::worldComm)
     && fileOperation::uniformFile(filePaths)
    );

    Pstream::broadcast(uniform, UPstream::worldComm);

    if (uniform)
    {
        if (Pstream::master(UPstream::worldComm))
        {
            if (!fName.empty())
            {
                IFstream is(fName);

                if (is.good())
                {
                    // Regular header or from decomposed data
                    ok = decomposedBlockData::readHeader(io, is);
                }
            }
        }

        Pstream::broadcasts
        (
            UPstream::worldComm,
            ok,
            io.headerClassName(),
            io.note()
        );
    }
    else
    {
        if (Pstream::nProcs(comm_) != Pstream::nProcs(UPstream::worldComm))
        {
            // Assume if different nprocs the communicators are also
            // different. Re-gather file paths on local master
            filePaths.resize(UPstream::nProcs(comm_));
            filePaths[UPstream::myProcNo(comm_)] = fName;
            Pstream::gatherList(filePaths, UPstream::msgType(), comm_);
        }

        // Intermediate storage arrays (master only)
        boolList result;
        wordList headerClassName;
        stringList note;

        if (Pstream::master(comm_))
        {
            const label np = Pstream::nProcs(comm_);

            result.resize(np, false);
            headerClassName.resize(np);
            note.resize(np);

            forAll(filePaths, proci)
            {
                if (!filePaths[proci].empty())
                {
                    if (proci > 0 && filePaths[proci] == filePaths[proci-1])
                    {
                        result[proci] = result[proci-1];
                        headerClassName[proci] = headerClassName[proci-1];
                        note[proci] = note[proci-1];
                    }
                    else
                    {
                        IFstream is(filePaths[proci]);

                        if (is.good())
                        {
                            result[proci] =
                                decomposedBlockData::readHeader(io, is);
                            headerClassName[proci] = io.headerClassName();
                            note[proci] = io.note();
                        }
                    }
                }
            }
        }

        // Is a more efficient scatter possible?
        PstreamBuffers pBufs(comm_);

        if (Pstream::master(comm_))
        {
            ok = result[0];
            io.headerClassName() = headerClassName[0];
            io.note() = note[0];

            // Scatter to each proc
            for (const int proci : pBufs.subProcs())
            {
                UOPstream os(proci, pBufs);
                os << result[proci] << headerClassName[proci] << note[proci];
            }
        }

        pBufs.finishedScatters();

        if (!Pstream::master(comm_))
        {
            UIPstream is(Pstream::masterNo(), pBufs);
            is >> ok >> io.headerClassName() >> io.note();
        }
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readHeader :" << " ok:" << ok
            << " class:" << io.headerClassName()
            << " for file:" << fName << endl;
    }
    return ok;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const word& typeName,
    const bool readOnProc
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readStream :"
            << " object : " << io.name()
            << " global : " << io.global()
            << " globalObject : " << io.globalObject()
            << " fName : " << fName << " readOnProc:" << readOnProc << endl;
    }

    // Close old stream
    io.close();

    autoPtr<ISstream> isPtr;
    bool isCollated = false;
    IOobject headerIO(io);

    // Detect collated format. This could be done on the local communicator
    // but we do it on the master node only for now.
    if (UPstream::master(UPstream::worldComm))
    {
        if (!fName.empty())
        {
            // This can happen in lagrangian field reading some processors
            // have no file to read from. This will only happen when using
            // normal writing since then the fName for the valid processors is
            // processorDDD/<instance>/.. . In case of collocated writing
            // the fName is already rewritten to processorsNN/.

            isPtr.reset(new IFstream(fName));

            if (isPtr->good())
            {
                // Read header data (on copy)
                headerIO.readHeader(*isPtr);

                isCollated = decomposedBlockData::isCollatedType(headerIO);

                if (!isCollated && !Pstream::parRun())
                {
                    // Short circuit: non-collated format. No parallel bits.
                    // Copy header and return.
                    if (debug)
                    {
                        Pout<< "masterUncollatedFileOperation::readStream :"
                            << " For object : " << io.name()
                            << " doing straight IFstream input from "
                            << fName << endl;
                    }

                    // Copy IOobject header information
                    io.IOobject::operator=(headerIO);
                    return isPtr;
                }
            }

            if (!isCollated)
            {
                // Close file. Reopened below.
                isPtr.clear();
            }
        }
    }

    Pstream::broadcast(isCollated, UPstream::worldComm);

    if (isCollated)
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream :"
                << " For object : " << io.name()
                << " starting collating input from " << fName << endl;
        }


        // Analyse the file path (on (co)master) to see the processors type
        // Note: this should really be part of filePath() which should return
        // both file and index in file.

        fileName path, procDir, local;
        procRangeType group;
        label nProcs;
        splitProcessorPath(fName, path, procDir, local, group, nProcs);


        if (!UPstream::parRun())
        {
            // Analyse the objectpath to find out the processor we're trying
            // to access
            label proci = detectProcessorPath(io.objectPath());

            if (proci == -1)
            {
                FatalIOErrorInFunction(*isPtr)
                    << "Could not detect processor number"
                    << " from objectPath:" << io.objectPath()
                    << " fileHandler : comm:" << comm_
                    << " ioRanks:" << flatOutput(ioRanks_)
                    << exit(FatalIOError);
            }

            // The local rank (offset)
            if (!group.empty())
            {
                proci = proci - group.start();
            }

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::readStream :"
                    << " For object : " << io.name()
                    << " starting input from block " << proci
                    << " of " << isPtr->name() << endl;
            }

            return decomposedBlockData::readBlock(proci, *isPtr, io);
        }
        else
        {
            // Are we reading from single-master file ('processors256') or
            // from multi-master files ('processors256_0-9')
            label readComm = -1;
            if (!group.empty())
            {
                readComm = comm_;
                if (UPstream::master(comm_) && !isPtr && !fName.empty())
                {
                    // In multi-master mode also open the file on the other
                    // masters
                    isPtr.reset(new IFstream(fName));

                    if (isPtr->good())
                    {
                        // Read header data (on copy)
                        IOobject headerIO(io);
                        headerIO.readHeader(*isPtr);
                    }
                }
            }
            else
            {
                // Single master so read on world
                readComm = UPstream::worldComm;
            }

            // Get size of file to determine communications type
            bool bigSize = false;

            if (Pstream::master(UPstream::worldComm))
            {
                // TBD: handle multiple masters?
                bigSize =
                (
                    off_t(Foam::fileSize(fName))
                  > off_t(maxMasterFileBufferSize)
                );
            }
            // Reduce (not broadcast)
            // - if we have multiple master files (FUTURE)
            Pstream::reduceOr(bigSize, UPstream::worldComm);

            const UPstream::commsTypes myCommsType
            (
                bigSize
              ? UPstream::commsTypes::scheduled
              : UPstream::commsTypes::nonBlocking
            );

            // Read my data
            return decomposedBlockData::readBlocks
            (
                readComm,
                fName,
                isPtr,
                io,
                myCommsType
            );
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream :"
                << " For object : " << io.name()
                << " starting separated input from " << fName << endl;
        }

        if (io.global() || io.globalObject())
        {
            // Use worldComm. Note: should not really need to gather filePaths
            // since we enforce sending from master anyway ...
            fileNameList filePaths(UPstream::nProcs(UPstream::worldComm));
            filePaths[UPstream::myProcNo(UPstream::worldComm)] = fName;
            Pstream::gatherList
            (
                filePaths,
                UPstream::msgType(),
                UPstream::worldComm
            );

            boolList readOnProcs
            (
                UPstream::listGatherValues<bool>
                (
                    readOnProc,
                    UPstream::worldComm
                )
            );
            // NB: local proc validity information required on sub-ranks too!
            readOnProcs.resize(UPstream::nProcs(UPstream::worldComm));
            readOnProcs[UPstream::myProcNo(UPstream::worldComm)] = readOnProc;

            // Uniform in local comm
            return read(io, UPstream::worldComm, true, filePaths, readOnProcs);
        }
        else
        {
            // Use local communicator
            fileNameList filePaths(UPstream::nProcs(comm_));
            filePaths[UPstream::myProcNo(comm_)] = fName;
            Pstream::gatherList
            (
                filePaths,
                UPstream::msgType(),
                comm_
            );

            boolList readOnProcs
            (
                UPstream::listGatherValues<bool>
                (
                    readOnProc,
                    comm_
                )
            );
            // NB: local proc validity information required on sub-ranks too!
            readOnProcs.resize(UPstream::nProcs(comm_));
            readOnProcs[UPstream::myProcNo(comm_)] = readOnProc;

            // Uniform in local comm
            const bool uniform = fileOperation::uniformFile(filePaths);

            return read(io, comm_, uniform, filePaths, readOnProcs);
        }
    }
}


bool Foam::fileOperations::masterUncollatedFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstreamOption::streamFormat format,
    const word& typeName
) const
{
    bool ok = true;

    if (io.global() || io.globalObject())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::read :"
                << " Reading global object " << io.name()
                << " worldComm:" << UPstream::worldComm
                << " Pstream::myProcNo:"
                << Pstream::myProcNo(UPstream::worldComm)
                << " amMaster:" << Pstream::master(UPstream::worldComm)
                << endl;
        }

        if (UPstream::master(UPstream::worldComm))
        {
            // Do master-only reading always.
            const bool oldParRun = UPstream::parRun(false);
            const int oldCache = fileOperation::cacheLevel(0);
            const label oldNProcs = fileOperation::nProcs();

            auto& is = io.readStream(typeName);
            ok = io.readData(is);
            io.close();

            const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
            fileOperation::cacheLevel(oldCache);
            UPstream::parRun(oldParRun);  // Restore parallel state
        }

        // Broadcast regIOobject content, with writeData/readData handling
        if (UPstream::parRun())
        {
            if (UPstream::master(UPstream::worldComm))
            {
                OPBstream os(UPstream::worldComm, format);

                os << io.headerClassName() << io.note();
                bool okWrite = io.writeData(os);
                ok = ok && okWrite;
            }
            else
            {
                IPBstream is(UPstream::worldComm, format);

                is >> io.headerClassName() >> io.note();
                ok = io.readData(is);
            }
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::read :"
                << " Reading local object " << io.name() << endl;
        }

        ok = io.readData(io.readStream(typeName));
        io.close();
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::read :"
            << " Read object:" << io.name()
            << " isGlobal:" << (io.global() || io.globalObject())
            << " status:" << ok << endl;
    }

    return ok;
}


bool Foam::fileOperations::masterUncollatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    fileName pathName(io.objectPath());

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::writeObject :"
            << " io:" << pathName << " writeOnProc:" << writeOnProc << endl;
    }

    // Make sure to pick up any new times
    setTime(io.time());

    // Update meta-data for current state
    const_cast<regIOobject&>(io).updateMetaData();

    autoPtr<OSstream> osPtr(NewOFstream(pathName, streamOpt, writeOnProc));
    OSstream& os = *osPtr;

    // If any of these fail, return (leave error handling to Ostream class)

    const bool ok =
    (
        os.good()
     && io.writeHeader(os)
     && io.writeData(os)
    );

    if (ok)
    {
        IOobject::writeEndDivider(os);
    }

    return ok;
}


Foam::instantList Foam::fileOperations::masterUncollatedFileOperation::findTimes
(
    const fileName& directory,
    const word& constantName
) const
{
    const auto iter = times_.cfind(directory);
    if (iter.good())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes :"
                << " Found " << iter.val()->size() << " cached times" << nl
                << "    for directory:" << directory << endl;
        }
        return *(iter.val());
    }
    else
    {
        instantList times;
        if (Pstream::master(UPstream::worldComm))
        {
            // Do master-only reading always.
            const bool oldParRun = UPstream::parRun(false);
            const int oldCache = fileOperation::cacheLevel(0);
            const label oldNProcs = fileOperation::nProcs();

            times = fileOperation::findTimes(directory, constantName);

            const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
            fileOperation::cacheLevel(oldCache);
            UPstream::parRun(oldParRun);  // Restore parallel state
        }

        Pstream::broadcast(times, UPstream::worldComm);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes :"
                << " Found times:" << flatOutput(times) << nl
                << "    for directory:" << directory << endl;
        }

        // Caching
        // - cache values even if no times were found since it might
        //   indicate a directory that is being filled later on ...
        if (cacheLevel() > 0)
        {
            auto* tPtr = new DynamicList<instant>(std::move(times));
            times_.set(directory, tPtr);

            return *tPtr;
        }

        // Times found (not cached)
        return times;
    }
}


void Foam::fileOperations::masterUncollatedFileOperation::setTime
(
    const Time& tm
) const
{
    if (tm.subCycling())
    {
        return;
    }

    // Mutable access to instant list for modification and sorting
    // - cannot use auto type deduction here

    auto iter = times_.find(tm.path());

    if (iter.good())
    {
        DynamicList<instant>& times = *(iter.val());

        const instant timeNow(tm.value(), tm.timeName());

        // The start index for checking and sorting (excluding "constant")
        const label startIdx =
        (
            (times.empty() || times[0].name() != tm.constant())
          ? 0
          : 1
        );

        // This routine always results in a sorted list of times, so first
        // check if the new time is greater than the latest existing time.
        // Can then simply append without extra searching or sorting

        if (times.size() <= startIdx || times.last() < timeNow)
        {
            times.append(timeNow);
        }
        else if
        (
            findSortedIndex
            (
                SubList<instant>(times, times.size()-startIdx, startIdx),
                timeNow
            ) < 0
        )
        {
            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::setTime :"
                    << " Caching time " << tm.timeName()
                    << " for case:" << tm.path() << endl;
            }

            times.append(timeNow);

            SubList<instant> realTimes
            (
                times, times.size()-startIdx, startIdx
            );
            Foam::stableSort(realTimes);
        }
    }

    fileOperation::setTime(tm);
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::NewIFstream
(
    const fileName& filePath
) const
{
    autoPtr<ISstream> isPtr;

    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the sub-ranks etc.

        fileNameList filePaths(Pstream::nProcs(comm_));
        filePaths[Pstream::myProcNo(comm_)] = filePath;
        Pstream::gatherList(filePaths, Pstream::msgType(), comm_);

        PstreamBuffers pBufs(comm_);

        if (Pstream::master(comm_))
        {
            // Same filename on the IO node -> same file
            const bool uniform = fileOperation::uniformFile(filePaths);

            if (uniform)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::NewIFstream :"
                        << " Opening global file " << filePath << endl;
                }

                readAndSend
                (
                    filePath,
                    identity(Pstream::nProcs(comm_)-1, 1),
                    pBufs
                );
            }
            else
            {
                for (const int proci : Pstream::subProcs(comm_))
                {
                    if (debug)
                    {
                        Pout<< "masterUncollatedFileOperation::NewIFstream :"
                            << " Opening local file " << filePath
                            << " for rank " << proci << endl;
                    }

                    readAndSend
                    (
                        filePaths[proci],
                        labelList(one{}, proci),
                        pBufs
                    );
                }
            }
        }


        pBufs.finishedSends();

        if (Pstream::master(comm_))
        {
            // Read myself
            isPtr.reset(new IFstream(filePaths[Pstream::masterNo()]));
        }
        else
        {
            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream :"
                    << " Reading " << filePath
                    << " from processor " << Pstream::masterNo() << endl;
            }

            List<char> buf(pBufs.recvDataCount(Pstream::masterNo()));

            if (!buf.empty())
            {
                UIPstream is(Pstream::masterNo(), pBufs);
                is.read(buf.data(), buf.size());
            }

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream :"
                    << " Done reading " << buf.size() << " bytes" << endl;
            }

            // A local character buffer copy of the Pstream contents.
            // Construct with same parameters (ASCII, current version)
            // as the IFstream so that it has the same characteristics.

            isPtr.reset(new ICharStream(std::move(buf)));

            // With the proper file name
            isPtr->name() = filePath;
        }
    }
    else
    {
        // Read myself
        isPtr.reset(new IFstream(filePath));
    }

    return isPtr;
}


Foam::autoPtr<Foam::OSstream>
Foam::fileOperations::masterUncollatedFileOperation::NewOFstream
(
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    return autoPtr<OSstream>
    (
        new masterOFstream
        (
            comm_,
            pathName,
            streamOpt,
            IOstreamOption::NO_APPEND,
            writeOnProc
        )
    );
}


Foam::autoPtr<Foam::OSstream>
Foam::fileOperations::masterUncollatedFileOperation::NewOFstream
(
    IOstreamOption::atomicType atomic,
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    return autoPtr<OSstream>
    (
        new masterOFstream
        (
            atomic,
            comm_,
            pathName,
            streamOpt,
            IOstreamOption::NO_APPEND,
            writeOnProc
        )
    );
}


void Foam::fileOperations::masterUncollatedFileOperation::flush() const
{
    fileOperation::flush();
    times_.clear();
}


void Foam::fileOperations::masterUncollatedFileOperation::sync()
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::sync :"
            << " syncing information across processors" << endl;
    }

    fileOperation::sync();


    wordList timeNames;
    List<DynamicList<instant>> instants;

    if (Pstream::master(UPstream::worldComm))
    {
        timeNames.resize(times_.size());
        instants.resize(times_.size());

        // Flatten into two lists to preserve key/val pairing
        label i = 0;
        forAllConstIters(times_, iter)
        {
            timeNames[i] = iter.key();
            instants[i] = std::move(*(iter.val()));
            ++i;
        }
    }

    Pstream::broadcasts(UPstream::worldComm, timeNames, instants);

    times_.clear();
    forAll(timeNames, i)
    {
        fileName dir(timeNames[i]);
        auto ptr = autoPtr<DynamicList<instant>>::New(std::move(instants[i]));

        if (Pstream::parRun() && !Pstream::master(UPstream::worldComm))
        {
            // Replace processor0 ending with processorDDD
            fileName path;
            fileName pDir;
            fileName local;
            procRangeType group;
            label numProcs;
            const label proci = splitProcessorPath
            (
                dir,
                path,
                pDir,
                local,
                group,
                numProcs
            );

            //Pout<< "**sync : From dir : " << dir << nl
            //    << "    path : " << path << nl
            //    << "    pDir : " << pDir << nl
            //    << "    local: " << local << nl
            //    << "    proci: " << proci << nl
            //    << endl;

            const label myProci = Pstream::myProcNo(UPstream::worldComm);

            if (proci != -1 && proci != myProci)
            {
                dir = path/"processor" + Foam::name(myProci);
            }
        }

        times_.insert(dir, ptr);
    }
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::addWatch
(
    const fileName& fName
) const
{
    label watchFd = -1;
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        watchFd = monitor().addWatch(fName);
    }

    Pstream::broadcast(watchFd, UPstream::worldComm);
    return watchFd;
}


bool Foam::fileOperations::masterUncollatedFileOperation::removeWatch
(
    const label watchIndex
) const
{
    bool ok = false;
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        ok = monitor().removeWatch(watchIndex);
    }

    Pstream::broadcast(ok, UPstream::worldComm);
    return ok;
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::findWatch
(
    const labelList& watchIndices,
    const fileName& fName
) const
{
    label index = -1;

    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        forAll(watchIndices, i)
        {
            if (monitor().getFile(watchIndices[i]) == fName)
            {
                index = i;
                break;
            }
        }
    }

    Pstream::broadcast(index, UPstream::worldComm);
    return index;
}


void Foam::fileOperations::masterUncollatedFileOperation::addWatches
(
    regIOobject& rio,
    const fileNameList& files
) const
{
    const labelList& watchIndices = rio.watchIndices();

    // Do on master and distribute effect to subprocs such that after
    // all have consistent numbering & files

    DynamicList<label> newWatchIndices;
    if (UPstream::master())
    {
        // Switch off comms inside findWatch/addWatch etc.
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);
        const label oldNProcs = fileOperation::nProcs();

        labelHashSet removedWatches(watchIndices);

        for (const fileName& f : files)
        {
            const label index = findWatch(watchIndices, f);

            if (index == -1)
            {
                newWatchIndices.push_back(addWatch(f));
            }
            else
            {
                // Existing watch
                newWatchIndices.push_back(watchIndices[index]);
                removedWatches.erase(index);
            }
        }

        // Remove any unused watches
        for (const label index : removedWatches)
        {
            removeWatch(watchIndices[index]);
        }

        const_cast<masterUncollatedFileOperation&>(*this).nProcs(oldNProcs);
        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);
    }
    Pstream::broadcast(newWatchIndices);

    rio.watchIndices() = newWatchIndices;
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::getFile
(
    const label watchIndex
) const
{
    fileName fName;
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        fName = monitor().getFile(watchIndex);
    }

    Pstream::broadcast(fName, UPstream::worldComm);
    return fName;
}


void Foam::fileOperations::masterUncollatedFileOperation::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        monitor().updateStates(true, false);
    }
}


Foam::fileMonitor::fileState
Foam::fileOperations::masterUncollatedFileOperation::getState
(
    const label watchFd
) const
{
    unsigned int state = fileMonitor::UNMODIFIED;
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        state = monitor().getState(watchFd);
    }

    Pstream::broadcast(state, UPstream::worldComm);
    return fileMonitor::fileState(state);
}


void Foam::fileOperations::masterUncollatedFileOperation::setUnmodified
(
    const label watchFd
) const
{
    if (!UPstream::parRun() || Pstream::master(UPstream::worldComm))
    {
        monitor().setUnmodified(watchFd);
    }
}


// ************************************************************************* //
