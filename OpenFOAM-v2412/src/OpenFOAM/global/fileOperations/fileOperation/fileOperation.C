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

#include "fileOperation.H"
#include "objectRegistry.H"
#include "labelIOList.H"
#include "registerSwitch.H"
#include "stringOps.H"
#include "Time.H"
#include "OSspecific.H"  // for Foam::isDir etc
#include <cinttypes>

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fileOperation, 0);
    defineRunTimeSelectionTable(fileOperation, word);
    defineRunTimeSelectionTable(fileOperation, comm);

    word fileOperation::defaultFileHandler
    (
        debug::optimisationSwitches().getOrAdd<word>
        (
            "fileHandler",
            //Foam::fileOperations::uncollatedFileOperation::typeName,
            "uncollated",
            keyType::LITERAL
        )
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum<Foam::fileOperation::pathType>
Foam::fileOperation::pathTypeNames_
({
    { fileOperation::NOTFOUND, "notFound" },
    { fileOperation::ABSOLUTE, "absolute" },
    { fileOperation::OBJECT, "objectPath" },
    { fileOperation::WRITEOBJECT, "writeObject" },
    { fileOperation::PROCUNCOLLATED, "uncollatedProc" },
    { fileOperation::PROCBASEOBJECT, "globalProc" },
    { fileOperation::PROCOBJECT, "localProc" },
    { fileOperation::PARENTOBJECT, "parentObjectPath" },
    { fileOperation::FINDINSTANCE, "findInstance" },
    { fileOperation::PROCUNCOLLATEDINSTANCE, "uncollatedProcInstance" },
    { fileOperation::PROCBASEINSTANCE, "globalProcInstance" },
    { fileOperation::PROCINSTANCE, "localProcInstance" }
});


Foam::word Foam::fileOperation::processorsBaseDir = "processors";

//- Caching (e.g. of time directories) - enabled by default
int Foam::fileOperation::cacheLevel_(1);

int Foam::fileOperation::nProcsFilter_(-1);

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Need to parse the numbers
// from "processors(\d+)" and
// from "processors(\d+)_(\d+)-(\d+)"
//
// Receive the string matching "^(\d+)(?:_(\d+)-(\d+))?/?$"
//
//    \1 = numProcs
//    \2 = firstProc
//    \3 = lastProc
//
// Return true on success and set parameters numProcs and group (size,start)
//
// Use low-level C-string to integer parsing to drive the sequence.
//
// For simplicity, also skip INT_MAX checks everywhere but check for
// - (errno) for success
// - (nptr == endptr) for leading junk
// - (*endptr != endChar) for trailing junk
// - skip INT_MAX checks as being too pessimistic

static bool parseProcsNumRange
(
    const std::string str,
    int& numProcs,
    Foam::fileOperation::procRangeType& group
)
{
    const char * nptr = str.c_str();
    char *endptr = nullptr;

    // 1. numProcs
    errno = 0;
    intmax_t parsed = std::strtoimax(nptr, &endptr, 10);
    if (errno || nptr == endptr) return false;  // bad parse

    const int nProcs = int(parsed);

    // End of string? Then no range and we are done.
    if (*endptr == '\0')
    {
        numProcs = nProcs;
        return true;
    }

    // Parse point at start of range ('_' character)?
    if (*endptr != '_') return false;
    nptr = ++endptr;


    // 2. firstProc
    errno = 0;
    parsed = std::strtoimax(nptr, &endptr, 10);
    if (errno || nptr == endptr) return false;  // bad parse

    const int firstProc = int(parsed);

    // Parse point at range separator ('-' character)?
    if (*endptr != '-') return false;
    nptr = ++endptr;


    // 3. lastProc
    errno = 0;
    parsed = std::strtoimax(nptr, &endptr, 10);
    if (errno || nptr == endptr) return false;  // bad parse

    const int lastProc = int(parsed);


    if
    (
        // Parse point at end of string
        (*endptr == '\0')

        // Input plausibility
        // Accept nProcs == 0 in case that becomes useful in the future
     && (nProcs >= 0 && firstProc >= 0 && firstProc <= lastProc)
    )
    {
        numProcs = nProcs;

        // Convert first/last to start/size
        group.reset(firstProc, lastProc-firstProc+1);

        return true;
    }

    return false;
}

} // End anonymous namespace


// Sorting of processor directories
#include "stringOpsSort.H"
namespace
{

// Sort processor directory names (natural order)
// - not strictly necessary
void sortProcessorDirs(Foam::UList<Foam::fileOperation::dirIndex>& dirs)
{
    if (dirs.size() > 1)
    {
        std::stable_sort
        (
            dirs.begin(),
            dirs.end(),
            []
            (
                const Foam::fileOperation::dirIndex& a,
                const Foam::fileOperation::dirIndex& b
            ) -> bool
            {
                return
                    Foam::stringOps::natural_sort::compare
                    (
                        a.first(),
                        b.first()
                    ) < 0;
            }
        );
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::instantList
Foam::fileOperation::sortTimes
(
    const fileNameList& dirEntries,
    const word& constantName
)
{
    // Check for "constant"
    bool haveConstant = false;

    if (!constantName.empty())
    {
        for (const fileName& dirName : dirEntries)
        {
            if (dirName == constantName)
            {
                haveConstant = true;
                break;
            }
        }
    }

    instantList times(dirEntries.size() + 1);
    label nTimes = 0;

    if (haveConstant)
    {
        times[nTimes].value() = 0;
        times[nTimes].name() = constantName;
        ++nTimes;
    }

    // Parse directory entries for scalar values
    for (const fileName& dirName : dirEntries)
    {
        if (readScalar(dirName, times[nTimes].value()))
        {
            times[nTimes].name() = dirName;
            ++nTimes;
        }
    }

    times.resize(nTimes);

    if (haveConstant)
    {
        if (nTimes > 2)
        {
            std::sort(&times[1], times.end(), instant::less());
        }
    }
    else if (nTimes > 1)
    {
        std::sort(times.begin(), times.end(), instant::less());
    }

    return times;
}


bool Foam::fileOperation::uniformFile(const fileNameList& names)
{
    if (names.empty())
    {
        return false;
    }

    const auto& object0 = names[0];

    for (label i = 1; i < names.size(); ++i)
    {
        if (object0 != names[i])
        {
            return false;
        }
    }

    return true;
}


bool Foam::fileOperation::uniformFile(const label comm, const fileName& name)
{
    if (!UPstream::parRun())
    {
        return true;
    }

    fileName masterName(name);

    Pstream::broadcast(masterName, comm);

    return returnReduceAnd((masterName == name), comm);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileMonitor& Foam::fileOperation::monitor() const
{
    if (!monitorPtr_)
    {
        monitorPtr_.reset
        (
            new fileMonitor
            (
                IOobject::fileModificationChecking == IOobject::inotify
             || IOobject::fileModificationChecking == IOobject::inotifyMaster
            )
        );
    }
    return *monitorPtr_;
}


void Foam::fileOperation::mergeTimes
(
    const instantList& extraTimes,
    const word& constantName,
    instantList& times
)
{
    if (extraTimes.size())
    {
        const bool haveConstant =
        (
            times.size()
         && times[0].name() == constantName
        );

        const bool haveExtraConstant =
        (
            extraTimes.size()
         && extraTimes[0].name() == constantName
        );

        // Combine times
        instantList combinedTimes(times.size()+extraTimes.size());
        label sz = 0;
        label extrai = 0;
        if (haveExtraConstant)
        {
            extrai = 1;
            if (!haveConstant)
            {
                combinedTimes[sz++] = extraTimes[0];    // constant
            }
        }
        forAll(times, i)
        {
            combinedTimes[sz++] = times[i];
        }
        for (; extrai < extraTimes.size(); extrai++)
        {
            combinedTimes[sz++] = extraTimes[extrai];
        }
        combinedTimes.setSize(sz);
        times.transfer(combinedTimes);

        // Sort
        if (times.size() > 1)
        {
            label starti = 0;
            if (times[0].name() == constantName)
            {
                starti = 1;
            }
            std::sort(&times[starti], times.end(), instant::less());

            // Filter out duplicates
            label newi = starti+1;
            for (label i = newi; i < times.size(); i++)
            {
                if (times[i].value() != times[i-1].value())
                {
                    if (newi != i)
                    {
                        times[newi] = times[i];
                    }
                    newi++;
                }
            }

            times.setSize(newi);
        }
    }
}


bool Foam::fileOperation::isFileOrDir(const bool isFile, const fileName& f)
{
    return (isFile ? Foam::isFile(f) : Foam::isDir(f));
}


Foam::refPtr<Foam::fileOperation::dirIndexList>
Foam::fileOperation::lookupAndCacheProcessorsPath
(
    const fileName& fName,
    const bool syncPar
) const
{
    // If path is local to a processor (e.g. contains 'processor2')
    // find the corresponding actual processor directory (e.g. 'processors4')
    // and index (2)

    // Behaviour affected by
    // - UPstream::parRun()
    // - syncPar : usually true, only uncollated does false. Determines
    //             if directory status gets synchronised
    // - distributed() : different processors have different roots
    // - fileModificationChecking : (uncollated only) do IO on master only
    // - nProcsFilter_ : if set to
    //       0 : accept any directory (e.g. for redistributePar where we don't
    //           know yet number of read/write procs)
    //      -1 : accept only processorsDDD where DDD is nProcs(worldComm)
    //      >0 : accept the exact mentioned number of prcessors


    // Collated : check whether/how to filter processorsXXX directory names
    const label targetNProcs
    (
        (UPstream::parRun() && nProcsFilter_ < 0)
      ? UPstream::nProcs(UPstream::worldComm)
      : nProcsFilter_
    );


    fileName path, pDir, local;
    procRangeType group;
    label numProcs;
    const label proci =
        splitProcessorPath(fName, path, pDir, local, group, numProcs);

    if (proci != -1)
    {
        const fileName procPath(path/pDir);

        if (cacheLevel() > 0)
        {
            const auto iter = procsDirs_.cfind(procPath);

            if (iter.good())
            {
                return iter.val();
            }
        }

        DynamicList<dirIndex> procDirs;
        fileNameList dirEntries;

        // Read all directories to see any beginning with processor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Note: use parallel synchronised reading so cache will be same
        //       order on all processors

        const bool readDirMasterOnly
        (
            UPstream::parRun() && !distributed()
         &&
            (
                IOobject::fileModificationChecking == IOobject::timeStampMaster
             || IOobject::fileModificationChecking == IOobject::inotifyMaster
            )
        );

        // The above selection excludes masterUncollated, which uses inotify or
        // timeStamp but provides its own internals for readDir() anyhow.

        if (readDirMasterOnly)
        {
            // Parallel and non-distributed
            // Read on master only and send to subProcs

            if (UPstream::master(UPstream::worldComm))
            {
                dirEntries = Foam::readDir(path, fileName::Type::DIRECTORY);

                DebugInfo
                    << "readDir on master: send " << dirEntries.size()
                    << " names to sub-processes" << endl;
            }

            Pstream::broadcast(dirEntries, UPstream::worldComm);
        }
        else
        {
            // Serial or distributed roots.
            // Handle readDir() with virtual method

            if (debug)
            {
                Pout<< "readDir without special master/send treatment"
                    << endl;
            }

            dirEntries = readDir(path, fileName::Type::DIRECTORY);
        }

        // Extract info from processorN or processorsNN
        // - highest processor number
        // - directory+offset containing data for proci

        label nProcs = 0;
        for (const fileName& dirN : dirEntries)
        {
            // Analyse directory name
            fileName rp, rd, rl;
            label rNum;
            const label readProci =
                splitProcessorPath(dirN, rp, rd, rl, group, rNum);

            nProcs = max(nProcs, readProci+1);

            Tuple2<pathType, int> pathTypeIdx(pathType::NOTFOUND, 0);

            if (proci == readProci)
            {
                // Found "processorN"
                pathTypeIdx.first() = pathType::PROCUNCOLLATED;
            }
            else if (rNum != -1)
            {
                if (targetNProcs > 1 && (targetNProcs != rNum))
                {
                    // Current directory can never contain wanted proc
                    //Pout<< "For fName:" << fName
                    //    << "Ignoring directory " << dirN
                    //    << " since parsed rNum:" << rNum
                    //    << " targetNProcs:" << targetNProcs
                    //   << endl;
                    //error::printStack(Pout);

                    continue;
                }

                // "processorsNN" or "processorsNN_start-end"
                nProcs = max(nProcs, rNum);

                if (group.empty())
                {
                    // "processorsNN"

                    if (proci < rNum || (nProcsFilter_ == 0))
                    {
                        // And it is also in range.
                        // Eg for "processors4": 3 is ok, 10 is not

                        pathTypeIdx.first() = pathType::PROCBASEOBJECT;
                        pathTypeIdx.second() = proci;
                    }
                }
                else if (group.contains(proci) || (nProcsFilter_ == 0))
                {
                    // "processorsNN_start-end"
                    // - save the local proc offset

                    pathTypeIdx.first() = pathType::PROCOBJECT;
                    pathTypeIdx.second() = (proci - group.start());
                }
            }

            if (pathTypeIdx.first() != pathType::NOTFOUND)
            {
                procDirs.append(dirIndex(dirN, pathTypeIdx));
            }
        }

        // Sort processor directory names (natural order)
        sortProcessorDirs(procDirs);


        // Global check of empty/exists.
        // 1 : empty directory
        // 2 : non-empty directory
        // 3 : mixed empty/non-empty directory (after reduce)
        // Combines andOp<bool>() and orOp<bool>() in single operation

        unsigned int procDirsStatus = (procDirs.empty() ? 1u : 2u);

        if (debug)
        {
            Pout<< "fileOperation::lookupProcessorsPath " << procPath
                << " detected:" << flatOutput(procDirs) << endl;
        }

        if (UPstream::parRun() && (!distributed() || syncPar))
        {
            reduce(procDirsStatus, bitOrOp<unsigned int>());  // worldComm

            if (procDirsStatus == 3u)
            {
                // Mixed empty/exists for procDirs.
                // Synthesize missing directory name (consistency in cache
                // existence).
                // Cannot reliably synthesize RANK-COLLATED, only COLLATED or
                // UNCOLLATED.
                //
                // RANK-COLLATED should have been read from its corresponding
                // master anyhow

                int flavour(pathType::PROCUNCOLLATED);
                for (const dirIndex& pDir : procDirs)
                {
                    flavour = max(flavour, int(pDir.second().first()));
                }

                reduce(nProcs, maxOp<label>());  // worldComm
                reduce(flavour, maxOp<int>());   // worldComm

                if (procDirs.empty())
                {
                    Tuple2<pathType, int> pathTypeIdx(pathType(flavour), 0);

                    if
                    (
                        pathTypeIdx.first() == pathType::PROCBASEOBJECT
                        // Do not restrict to currently used processors
                        // && proci < nProcs
                    )
                    {
                        pathTypeIdx.second() = proci;

                        procDirs.append
                        (
                            dirIndex
                            (
                                processorsBaseDir + Foam::name(nProcs),
                                pathTypeIdx
                            )
                        );
                    }
                    else
                    {
                        // - pathType::PROCUNCOLLATED
                        // - poor fallback for pathType::PROCOBJECT
                        // - out-of-range pathType::PROCBASEOBJECT

                        procDirs.append
                        (
                            dirIndex
                            (
                                "processor" + Foam::name(proci),
                                pathTypeIdx
                            )
                        );
                    }

                    if (debug)
                    {
                        Pout<< "fileOperation::lookupProcessorsPath "
                            << procPath
                            << " synthetic:" << procDirs << endl;
                    }
                }
            }
        }
        else if (!UPstream::parRun())
        {
            // Serial: use the number of decompositions (if found)
            if (nProcs)
            {
                const_cast<fileOperation&>(*this).nProcs(nProcs);
            }
        }

        if (procDirsStatus & 2u)
        {
            if (cacheLevel() > 0)
            {
                procsDirs_.insert(procPath, procDirs);

                // Make sure to return a reference
                return procsDirs_[procPath];
            }
            else
            {
                return refPtr<dirIndexList>::New(procDirs);
            }
        }
    }

    return refPtr<dirIndexList>::New();
}


Foam::refPtr<Foam::fileOperation::dirIndexList>
Foam::fileOperation::lookupProcessorsPath(const fileName& fName) const
{
    return lookupAndCacheProcessorsPath(fName, true);
}


bool Foam::fileOperation::exists(IOobject& io) const
{
    // Generate output filename for object
    fileName objPath(objectPath(io, word::null));

    // Test for either directory or a (valid) file & IOobject
    bool ok;
    if (io.name().empty())
    {
        ok = isDir(objPath);
    }
    else
    {
        ok =
        (
            isFile(objPath)
            // object with local scope
         && io.typeHeaderOk<labelIOList>(false)
        );
    }

    if (!ok)
    {
        // Re-test with searched for objectPath. This is for backwards
        // compatibility
        fileName originalPath(filePath(io.objectPath()));
        if (originalPath != objPath)
        {
            // Test for either directory or a (valid) file & IOobject
            if (io.name().empty())
            {
                ok = isDir(originalPath);
            }
            else
            {
                ok =
                (
                    isFile(originalPath)
                    // object with local scope
                 && io.typeHeaderOk<labelIOList>(false)
                );
            }
        }
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperation::fileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots
)
:
    comm_(commAndIORanks.first()),
    nProcs_(UPstream::nProcs(UPstream::worldComm)),
    distributed_(distributedRoots),
    ioRanks_(commAndIORanks.second())
{}


Foam::fileOperation::fileOperation
(
    const label comm,
    const labelUList& ioRanks,
    const bool distributedRoots
)
:
    comm_(comm),
    nProcs_(UPstream::nProcs(UPstream::worldComm)),
    distributed_(distributedRoots),
    ioRanks_(ioRanks)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    return io.objectPath();
}


bool Foam::fileOperation::writeObject
(
    const regIOobject& io,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (writeOnProc)
    {
        const fileName pathName(io.objectPath());

        mkDir(pathName.path());

        autoPtr<OSstream> osPtr(NewOFstream(pathName, streamOpt));

        if (!osPtr)
        {
            return false;
        }

        OSstream& os = *osPtr;

        // Update meta-data for current state
        const_cast<regIOobject&>(io).updateMetaData();

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
    return true;
}


Foam::fileName Foam::fileOperation::filePath
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    if (debug)
    {
        Pout<< "fileOperation::filePath :" << " fName:" << fName << endl;
    }

    fileName path, pDir, local;
    procRangeType group;
    label numProcs;
    label proci =
        splitProcessorPath(fName, path, pDir, local, group, numProcs);

    if (numProcs != -1)
    {
        WarningInFunction << "Filename is already adapted:" << fName << endl;
    }

    // Give preference to processors variant
    fileName foundName;
    if (proci != -1)
    {
        // Get all processor directories
        refPtr<dirIndexList> procDirs(lookupProcessorsPath(fName));
        for (const dirIndex& dirIdx : procDirs())
        {
            const fileName& procDir = dirIdx.first();

            fileName collatedName(path/procDir/local);
            if (exists(collatedName, checkGzip, followLink))
            {
                if (debug)
                {
                    Pout<< "fileOperation::filePath : " << collatedName << endl;
                }
                foundName = collatedName;
            }
        }
    }

    //if (returnReduceOr(foundName.empty()))    // worldComm
    if (foundName.empty())
    {
        // There is at least one processor that cannot find the processors
        // directory. Re-do with straight supplied filename
        if (exists(fName, checkGzip, followLink))
        {
            if (foundName.empty())
            {
                foundName = fName;
            }
        }
    }

    if (!foundName.empty())
    {
        if (debug)
        {
            Pout<< "fileOperation::filePath : " << foundName << endl;
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "fileOperation::filePath : Not found" << endl;
        }
    }
    return foundName;
}


Foam::label Foam::fileOperation::addWatch(const fileName& fName) const
{
    return monitor().addWatch(fName);
}


bool Foam::fileOperation::removeWatch(const label watchIndex) const
{
    return monitor().removeWatch(watchIndex);
}


Foam::label Foam::fileOperation::findWatch
(
    const labelList& watchIndices,
    const fileName& fName
) const
{
    forAll(watchIndices, i)
    {
        if (getFile(watchIndices[i]) == fName)
        {
            return i;
        }
    }
    return -1;
}


void Foam::fileOperation::addWatches
(
    regIOobject& rio,
    const fileNameList& files
) const
{
    const labelList& watchIndices = rio.watchIndices();

    DynamicList<label> newWatchIndices;
    labelHashSet removedWatches(watchIndices);

    for (const fileName& f : files)
    {
        const label index = findWatch(watchIndices, f);

        if (index == -1)
        {
            newWatchIndices.append(addWatch(f));
        }
        else
        {
            // Existing watch
            newWatchIndices.append(watchIndices[index]);
            removedWatches.erase(index);
        }
    }

    // Remove any unused watches
    for (const label index : removedWatches)
    {
        removeWatch(watchIndices[index]);
    }

    rio.watchIndices() = newWatchIndices;
}


Foam::fileName Foam::fileOperation::getFile(const label watchIndex) const
{
    return monitor().getFile(watchIndex);
}


void Foam::fileOperation::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    monitor().updateStates(masterOnly, UPstream::parRun());
}


Foam::fileMonitor::fileState Foam::fileOperation::getState
(
    const label watchFd
) const
{
    return monitor().getState(watchFd);
}


void Foam::fileOperation::setUnmodified(const label watchFd) const
{
    monitor().setUnmodified(watchFd);
}


Foam::instantList Foam::fileOperation::findTimes
(
    const fileName& directory,
    const word& constantName
) const
{
    if (debug)
    {
        Pout<< "fileOperation::findTimes : Finding times in directory "
            << directory << endl;
    }

    // Note: do NOT use master-only reading here (as per lookupProcessorsPath)
    // since this routine is called on an individual processorN directory

    // Read directory entries into a list
    fileNameList dirEntries(Foam::readDir(directory, fileName::DIRECTORY));
    instantList times = sortTimes(dirEntries, constantName);


    // Get all processor directories
    refPtr<dirIndexList> procDirs(lookupProcessorsPath(directory));
    for (const dirIndex& dirIdx : procDirs())
    {
        const fileName& procDir = dirIdx.first();
        fileName collDir(processorsPath(directory, procDir));
        if (!collDir.empty() && collDir != directory)
        {
            fileNameList extraEntries
            (
                Foam::readDir
                (
                    collDir,
                    fileName::DIRECTORY
                )
            );
            mergeTimes
            (
                sortTimes(extraEntries, constantName),
                constantName,
                times
            );
        }
    }

    if (debug)
    {
        Pout<< "fileOperation::findTimes : Found times:" << flatOutput(times)
            << endl;
    }
    return times;
}


Foam::IOobject Foam::fileOperation::findInstance
(
    const IOobject& startIO,
    const scalar startValue,
    const word& stopInstance,
    const bool constant_fallback
) const
{
    const Time& time = startIO.time();
    IOobject io(startIO);

    // Note: - if name is empty, just check the directory itself
    //       - check both for isFile and headerOk since the latter does a
    //         filePath so searches for the file.
    //       - check for an object with local file scope (so no looking up in
    //         parent directory in case of parallel)

    if (exists(io))
    {
        DebugInFunction
            << "Found exact match for \"" << io.name()
            << "\" in " << io.instance()/io.local()
            << endl;

        return io;
    }

    // Handling failures afterwards
    const bool exitIfMissing = startIO.isReadRequired();

    enum failureCodes { FAILED_STOPINST = 1, FAILED_CONSTINST = 2 };
    int failed(0);

    instantList ts = time.times();

    {
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
            io.instance() = ts[instIndex].name();

            // Shortcut: if actual directory is the timeName we've
            // already tested it
            if
            (
                io.instance() == startIO.instance()
             && io.instance() != stopInstance
            )
            {
                continue;
            }

            if (exists(io))
            {
                DebugInFunction
                    << "Found exact match for \"" << io.name()
                    << "\" in " << io.instance()/io.local()
                    << endl;

                return io;
            }

            // Check if hit minimum instance
            if (io.instance() == stopInstance)
            {
                DebugInFunction
                    << "Hit stopInstance " << stopInstance << endl;

                if (exitIfMissing)
                {
                    failed = failureCodes::FAILED_STOPINST;
                }
                else
                {
                    // At the stopInstance
                    return io;
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
            !failed
         && (ts.empty() || ts[0].name() != time.constant() || startValue < 0)
        )
        {
            io.instance() = time.constant();

            if (exists(io))
            {
                DebugInFunction
                    << "Found constant match for \"" << io.name()
                    << "\" in " << io.instance()/io.local()
                    << endl;

                return io;
            }
        }
    }


    if (!failed)
    {
        if (exitIfMissing)
        {
            failed = failureCodes::FAILED_CONSTINST;
        }
        else if (constant_fallback)
        {
            io.instance() = time.constant();
        }
        else
        {
            io.instance().clear();
        }
    }


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

    return io;
}


Foam::fileNameList Foam::fileOperation::readObjects
(
    const objectRegistry& obr,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "fileOperation::readObjects :"
            << " object-path:" << obr.objectPath()
            << " instance" << instance
            << " local:" << local << endl;
    }

    // dbDir() is relative to Time,
    // so use dbDir() from self or from parent, but not both!
    fileName path;
    if (obr.dbDir().empty() || !obr.parent().dbDir().empty())
    {
        path = obr.path(instance, local);
    }
    else
    {
        path = obr.path(instance, obr.dbDir()/local);
    }

    newInstance.clear();
    fileNameList objectNames;

    if (Foam::isDir(path))
    {
        newInstance = instance;
        objectNames = Foam::readDir(path, fileName::Type::FILE);
    }
    else
    {
        // Get processors equivalent of path
        fileName procsPath(filePath(path));

        if (!procsPath.empty())
        {
            newInstance = instance;
            objectNames = Foam::readDir(procsPath, fileName::Type::FILE);
        }
    }
    return objectNames;
}


Foam::label Foam::fileOperation::nProcs
(
    const fileName& dir,
    const fileName& local,
    const label wanted      // expected nProcs. 0 if not supplied
) const
{
    label nProcs = 0;
    if (UPstream::master(comm_))
    {
        fileNameList dirNames(Foam::readDir(dir, fileName::Type::DIRECTORY));

        // E.g. (maybe left over)
        // processor0
        // processor1
        // processors5
        // processor90     // gap!
        // processors10
        //
        //- if wanted is 2 then processor0,1 is enough
        //- if wanted is not set then return highest valid range
        //  (processors10 in above case)
        //- if wanted cannot be matched (e.g. 37) return 0

        // For marking off contiguous ranges
        bitSet foundDirs;
        if (wanted > 0)
        {
            foundDirs.resize(wanted);
        }

        // Detect any processorsDDD or processorDDD
        label maxProc = -1;
        for (const fileName& dirN : dirNames)
        {
            fileName rp, rd, rl;
            procRangeType group;
            label rNum;

            const label readProci =
                splitProcessorPath(dirN, rp, rd, rl, group, rNum);

            maxProc = max(maxProc, readProci);

            if (rNum > 0)       // processorsDDD where DDD>0
            {
                // Use processors number
                maxProc = max(maxProc, rNum-1);
                // Mark in cache. TBD separate handling for groups?
                foundDirs.set(labelRange(0, rNum));

                if (wanted == rNum)
                {
                    // Exact match of processorsDDD. Direct exit.
                    maxProc = rNum-1;
                    foundDirs.resize(rNum);
                    break;
                }
            }
            else if (readProci >= 0)
            {
                // Mark in cache
                foundDirs.set(readProci);
            }
        }
        nProcs = maxProc+1;

        // Override with any gaps in processorDDD numbering (can never happen
        // with collated)
        const label gapIndex = foundDirs.find_first_not();
        if (gapIndex > 0)
        {
            nProcs = gapIndex-1;
        }

        if (nProcs == 0 && Foam::isDir(dir/processorsBaseDir))
        {
            WarningInFunction
                << "Defunct collated naming: " << processorsBaseDir << nl
                << "Manually rename with the decomposition number. Eg,"
                << nl << nl
                << "    mv processors processors16" << nl << nl
                << "...returning 1" << endl;

            nProcs = 1;
        }
    }
    Pstream::broadcast(nProcs, comm_);
    return nProcs;
}


void Foam::fileOperation::flush() const
{
    if (debug)
    {
        Pout<< "fileOperation::flush : clearing processor directories cache"
            << endl;
    }
    procsDirs_.clear();
}


void Foam::fileOperation::sync()
{
    if (debug)
    {
        Pout<< "fileOperation::sync : parallel synchronisation"
            << endl;
    }

    Pstream::broadcasts
    (
        UPstream::worldComm,
        nProcs_,
        procsDirs_
    );
}


Foam::fileName Foam::fileOperation::processorsCasePath
(
    const IOobject& io,
    const word& procsDir
) const
{
    return io.rootPath()/io.globalCaseName()/procsDir;
}


Foam::fileName Foam::fileOperation::processorsPath
(
    const IOobject& io,
    const word& instance,
    const word& procsDir
) const
{
    return
        processorsCasePath(io, procsDir)
       /instance
       /io.db().dbDir()
       /io.local();
}


Foam::fileName Foam::fileOperation::processorsPath
(
    const fileName& dir,
    const word& procsDir
) const
{
    // Check if directory is processorDDD

    const word caseName(dir.name());
    if (caseName.starts_with("processor"))
    {
        // Reject both '^processor$' and '^processors.*$'

        if (!std::isdigit(caseName[9]))
        {
            WarningInFunction << "Directory " << dir
                << " does not end in old-style processorDDD" << endl;
        }

        return dir.path()/procsDir;
    }

    return fileName();
}


Foam::label Foam::fileOperation::splitProcessorPath
(
    const fileName& objPath,
    fileName& path,
    fileName& procDir,
    fileName& local,

    procRangeType& group,
    label& nProcs
)
{
    // Return value
    label returnProci = -1;

    // Clear out the return parameters

    path.clear();
    procDir.clear();
    local.clear();
    group.clear();

    // Invalidate detected number of processors
    nProcs = -1;

    // The local processor group is read as first/last, but stored as
    // start/size.  Empty with start=0, size=0 if no range is detected


    // Start of 'processor..' directory name (the procDir)
    size_t pos = 0;

    // The slash starting the trailing (local) directory
    size_t slashLocal = string::npos;


    // Search for processor at start of string or after /processor
    //
    // 'processor(\d+)'
    // 'processors(\d+)'
    // 'processors(\d+)_(\d+)-(\d+)'

    for
    (
        /*nil*/;
        (pos = objPath.find("processor", pos)) != string::npos;
        pos += 9
    )
    {
        if (pos > 0 && objPath[pos-1] != '/')
        {
            // Not start of string or after /processor
            continue;
        }

        // The parse point. One past 'processor'
        size_t firstp = pos + 9;

        // normal: 'processor(\d+)'
        // plural: 'processors(\d+)'

        const bool plural = (objPath[firstp] == 's');

        if (plural)
        {
            ++firstp;  // Skip over the 's'
        }
        else if (!std::isdigit(objPath[firstp]))
        {
            // Non-plural version (uncollated) requires digits only
            continue;
        }

        // The next slash indicates there is a local directory
        slashLocal = objPath.find('/', firstp);

        // The last parse point is the slash, or end of string
        const size_t lastp =
            (slashLocal == string::npos ? objPath.length() : slashLocal);

        if (!std::isdigit(objPath[lastp-1]))
        {
            // Must end in a digit!
            // This traps entries that are too short or look quite wrong
            // and avoid a string to int conversion that will fail anyhow
            continue;
        }


        // Match: '^processors(\d+)$'  -> nProcs

        // Match: '^processors(\d+)_(\d+)-(\d+)$'
        // \1 = nProcs
        // \2 = beg processor group
        // \3 = end processor group (inclusive)

        if (plural)
        {
            int nProcsRead = 0;

            if
            (
                parseProcsNumRange
                (
                    objPath.substr(firstp, lastp-firstp),
                    nProcsRead,
                    group
                )
            )
            {
                // Total number of processors
                nProcs = nProcsRead;

                // We are done!
                break;
            }
        }

        // Single
        // Match: '^processor(\d+)$'   -> proci

        label proci = 0;
        if
        (
            Foam::read(objPath.substr(firstp, lastp-firstp), proci)
         && (proci >= 0)
        )
        {
            // Capture value of an individual processor
            returnProci = proci;

            // We are done!
            break;
        }
    }

    if (pos != string::npos)
    {
        // The split succeeded, extract the components.

        // The leading directory
        if (pos > 0)
        {
            path = objPath.substr(0, pos-1);
        }

        // The slash starting the trailing (local) directory
        if (slashLocal != string::npos)
        {
            procDir = objPath.substr(pos, slashLocal-pos);
            local = objPath.substr(slashLocal+1);
        }
        else
        {
            procDir = objPath.substr(pos);
        }
    }

    return returnProci;
}


Foam::label Foam::fileOperation::detectProcessorPath(const fileName& fName)
{
    fileName path, pDir, local;
    procRangeType group;
    label nProcs;
    return splitProcessorPath(fName, path, pDir, local, group, nProcs);
}


// * * * * * * * * * * * * * *  Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<fileOperation>& iproxy
)
{
    const auto& fp = *iproxy;

    os  << "fileHandler:" << fp.type()
        << " nProcs:" << fp.nProcs()
        << " comm:" << fp.comm()
        << " distributed:" << fp.distributed()
        << " ioranks: " << flatOutput(fp.ioRanks())
        << " ranks: ";

    if (fp.comm() >= 0)
    {
        os << flatOutput(UPstream::procID(fp.comm()));
    }
    os  << nl;

    return os;
}

// ************************************************************************* //
