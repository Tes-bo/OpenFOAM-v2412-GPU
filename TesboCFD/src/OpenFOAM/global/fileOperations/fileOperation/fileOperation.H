/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

Class
    Foam::fileOperation

Description
    An encapsulation of filesystem-related operations.

    Several of the file handlers can be configured to use specific ranks for
    IO operations. These can either be defined from the command-line
    \c -ioRanks option, or via the \c FOAM_ENV environment variable.
    In either case, the list of IO ranks shall always include the value
    \b 0 (master rank). The ranks may be specified as a plain list, or as an
    OpenFOAM list. The special value \em host can be used to specify selection
    based on hostname.

    For example,
    \verbatim
        FOAM_IORANKS='0 4 8'   decomposePar -fileHandler hostCollated
        FOAM_IORANKS='0,4,8'   decomposePar -fileHandler hostCollated
        FOAM_IORANKS='(0 4 8)' decomposePar -fileHandler hostCollated
    \endverbatim
    will generate
    \verbatim
        processors12_0-3/   : containing data for processors 0 to 3
        processors12_4-7/   : containing data for processors 4 to 7
        processors12_8-11/  : containing data for processors 8 to 11
    \endverbatim

    The value \em host can be used to specify a single IO rank for each
    host. For example,
    \verbatim
        decomposePar -fileHandler collated -ioRanks host
        decomposePar -fileHandler hostCollated
    \endverbatim

Environment
    - \c FOAM_ENV : list of io-ranks as plain space or comma separated
      list or as an OpenFOAM formatted list. Eg, '(0 4 8)'

Namespace
    Foam::fileOperations

Description
    Namespace for implementations of a fileOperation

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperation_H
#define Foam_fileOperation_H

#include "ISstream.H"
#include "Ostream.H"
#include "UPstream.H"
#include "fileMonitor.H"
#include "fileNameList.H"
#include "instantList.H"
#include "refPtr.H"
#include "bitSet.H"
#include "Enum.H"
#include "Tuple2.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fileOperation;
class objectRegistry;
class regIOobject;
class IOobject;
class Time;

Ostream& operator<<(Ostream& os, const InfoProxy<fileOperation>& info);

/*---------------------------------------------------------------------------*\
                         Class fileOperation Declaration
\*---------------------------------------------------------------------------*/

class fileOperation
{
public:

    // Public Data Types

        //- Enumeration for the location of an IOobject
        enum pathType : int
        {
            NOTFOUND = 0,           //!< Not found
            ABSOLUTE,               //!< instance is absolute directory
            OBJECT,                 //!< io.objectPath() exists
            WRITEOBJECT,            //!< write path exists

            // NOTE: increasing precedence (uncollated, collated, rank-collated)

            PROCUNCOLLATED,
                        //!< objectPath exists in 'processorN'
            PROCBASEOBJECT = PROCUNCOLLATED + 1,
                        //!< objectPath exists in 'processorsNN'
            PROCOBJECT = PROCBASEOBJECT + 1,
                        //!< objectPath exists in 'processorsNN_first-last'

            PARENTOBJECT,           //!< parent of object path
            FINDINSTANCE,           //!< file found in time directory
            PROCUNCOLLATEDINSTANCE, //!< as PROCUNCOLLATED but with instance
            PROCBASEINSTANCE,       //!< as PROCBASEOBJECT but with instance
            PROCINSTANCE            //!< as PROCOBJECT but with instance
        };
        static const Enum<pathType> pathTypeNames_;

        //- Augment fileName with pathType and local offset
        typedef Tuple2<fileName, Tuple2<pathType, int>> dirIndex;
        typedef List<dirIndex> dirIndexList;

        //- For addressing a range of processors,
        //- identical to UPstream::rangeType
        typedef IntRange<int> procRangeType;


private:

    // Private Static Data

        //- Storage of the dummy file handler (demand-driven)
        static refPtr<fileOperation> dummyHandlerPtr_;

        //- Filtering of processorsDDD directories (in parallel):
        //      - >0 : only accept processorsDDD (default)
        //      -  0 : accept anything (e.g. when detecting nprocs from
        //             directory naming)
        //      - -1 : use nProcs(UPstream::worldComm)
        static int nProcsFilter_;


protected:

    // Protected Static Data

        //- Cache level (eg, for caching time directories). Default: 1
        static int cacheLevel_;


    // Protected Data

        //- Communicator to use
        mutable label comm_;

        //- Overall number of processors.
        //  Used to synthesise processor directory naming:
        //  - parallel: UPstream::nProcs(UPstream::commWorld())
        //  - non-parallel: detected from processor dir naming ('processorsNN')
        label nProcs_;

        //- Distributed roots (parallel run)
        mutable bool distributed_;

        //- The list of IO ranks (global ranks)
        //  Primarily for additional bookkeeping in non-parallel
        const labelList ioRanks_;

        //- Detected processors directories
        mutable HashTable<dirIndexList> procsDirs_;

        //- File-change monitor for all registered files
        mutable std::unique_ptr<fileMonitor> monitorPtr_;


   // Protected Member Functions

        //- Get or create fileMonitor singleton
        fileMonitor& monitor() const;

        //- Merge two times
        static void mergeTimes
        (
            const instantList& extraTimes,
            const word& constantName,
            instantList& times
        );

        //- Helper: check for file (isFile) or directory (!isFile)
        static bool isFileOrDir(const bool isFile, const fileName&);

        //- Lookup name of processorsDDD using cache.
        //  \return empty fileName if not found.
        refPtr<dirIndexList> lookupAndCacheProcessorsPath
        (
            const fileName& objectPath,
            const bool syncPar
        ) const;

        //- Lookup name of processorsDDD using cache.
        //  \note To be called on all processors
        //  \return empty fileName if not found.
        virtual refPtr<dirIndexList> lookupProcessorsPath
        (
            const fileName& objectPath
        ) const;

        //- Does IOobject exist?
        //- Is either a directory (empty name()) or a file
        bool exists(IOobject& io) const;


        //- Is proci a master rank in the communicator (in parallel)
        //- or a master rank in the IO ranks (non-parallel)
        bool isIOrank(const label proci) const;

        //- Helper: output which ranks are IO
        void printRanks() const;


        //- Construction helper: check for locally allocated communicator
        static inline label getManagedComm(const label communicator)
        {
            return
            (
                (
                    communicator < 0
                 || communicator == UPstream::commGlobal()
                 || communicator == UPstream::commSelf()
                 || communicator == UPstream::commWorld()
                )
              ? -1
              : communicator
            );
        }


private:

    // Private Member Functions

        //- Low-level implementation for subsetting 'clone' operation
        static autoPtr<fileOperation> New_impl
        (
            const fileOperation& origHandler,
            //! The selected sub-ranks (in comm-world)
            const labelUList& subProcs,
            bool verbose
        );


public:

    // Static Data

        //- Return the processors directory name (usually "processors")
        static word processorsBaseDir;

        //- Name of the default fileHandler
        static word defaultFileHandler;


    // Public Data Types

        //- Runtime type information
        TypeName("fileOperation");


    //- The currently active file handler. Avoid accessing directly
    static refPtr<fileOperation> fileHandlerPtr_;

    //- Reference to a dummy file handler.
    static refPtr<fileOperation> null();


    // Constructors

        //- Construct from communicator,
        //- optionally with specified io-ranks and/or distributed roots
        explicit fileOperation
        (
            const label comm,
            const labelUList& ioRanks = labelUList::null(),
            const bool distributedRoots = false
        );

        //- Construct from communicator with specified io-ranks
        explicit fileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots = false
        );


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperation,
            word,
            (
                bool verbose
            ),
            (verbose)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperation,
            comm,
            (
                const Tuple2<label, labelList>& commAndIORanks,
                const bool distributedRoots,
                bool verbose
            ),
            (commAndIORanks, distributedRoots, verbose)
        );


    // Selectors

        //- Select fileHandler-type.
        //- Uses defaultFileHandler if the handlerType is empty.
        static autoPtr<fileOperation> New
        (
            const word& handlerType,
            bool verbose = false
        );

        //- Select fileHandler-type.
        //- Uses defaultFileHandler if the handlerType is empty.
        static autoPtr<fileOperation> New
        (
            const word& handlerType,
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


        //- Replicate the given fileHandler properties with
        //- a subset of (global) ranks.
        //  Always includes rank 0 and constrained by the global numProcs
        //  Returns nullptr if myProcNo is not involved.
        static autoPtr<fileOperation> New
        (
            //! The original file handler
            const fileOperation& origHandler,
            //! The selected sub-ranks (in comm-world)
            const boolUList& useProc,
            bool verbose = false
        );

        //- Replicate the given fileHandler properties with
        //- a subset of (global) ranks.
        //  Always includes rank 0 and constrained by the global numProcs
        //  Returns nullptr if myProcNo is not involved.
        static autoPtr<fileOperation> New
        (
            //! The original file handler
            const fileOperation& origHandler,
            //! The selected sub-ranks (in comm-world). Monotonic order!
            const bitSet& useProc,
            bool verbose = false
        );


    //- Destructor
    virtual ~fileOperation() = default;


   // Factory Methods, Singleton-type Functions

        //- The commonly used uncollatedFileOperation
        static autoPtr<fileOperation> NewUncollated();

        //- Return the current file handler.
        //- Will create the default file handler if necessary.
        static const fileOperation& fileHandler();

        //- Delete current file handler.
        //  \returns the old handler.
        //  Should have [[nodiscard]], but gcc ignores void casting.
        static refPtr<fileOperation> fileHandler(std::nullptr_t);

        //- Replace the current file handler.
        //  The following are considered no-ops:
        //  - an empty/invalid newHandler does \b not delete, use a literal
        //    \c nullptr (std::nullptr_t) for that
        //  - if new handler and current handler are identical (same pointer).
        //  .
        //  \returns the old handler (on change), nullptr otherwise
        //  Should have [[nodiscard]], but gcc ignores void casting.
        static refPtr<fileOperation> fileHandler
        (
            refPtr<fileOperation>& newHandler
        );

        //- Replace the current file handler.
        //  The following are considered no-ops:
        //  - an empty/invalid newHandler does \b not delete, use a literal
        //    \c nullptr (std::nullptr_t) for that
        //  - if new handler and current handler are identical (same pointer).
        //  .
        //  \returns the old handler (on change), nullptr otherwise
        //  Should have [[nodiscard]], but gcc ignores void casting.
        static refPtr<fileOperation> fileHandler
        (
            refPtr<fileOperation>&& newHandler
        );

        //- Replace the current file handler.
        //  The following are considered no-ops:
        //  - an empty/invalid newHandler does \b not delete, use a literal
        //    \c nullptr (std::nullptr_t) for that
        //  - if new handler and current handler are identical (same pointer).
        //  .
        //  \returns the old handler (on change), nullptr otherwise
        //  Should have [[nodiscard]], but gcc ignores void casting.
        static refPtr<fileOperation> fileHandler
        (
            autoPtr<fileOperation>&& newHandler
        );


   // Static Functions

        //- Return cache level
        static int cacheLevel() noexcept
        {
            return cacheLevel_;
        }

        //- Set cache level (0 = off). \return the previous value
        static int cacheLevel(int level) noexcept
        {
            int old(cacheLevel_);
            cacheLevel_ = level;
            return old;
        }

        //- Return collated 'processorsDDD' filtering
        static int nProcsFilter() noexcept
        {
            return nProcsFilter_;
        }

        //- Set collated 'processorsDDD' filtering (0 = off).
        //  \return the previous value
        static int nProcsFilter(int level) noexcept
        {
            int old(nProcsFilter_);
            nProcsFilter_ = level;
            return old;
        }

        //- Sort directory entries according to time value,
        //  with "constant" appearing first (if it exists)
        static instantList sortTimes
        (
            const fileNameList& dirEntries,
            const word& constantName = "constant"
        );

        //- True if the file names are identical. False on an empty list
        static bool uniformFile(const fileNameList& names);

        //- True if the file name is identical on all ranks
        static bool uniformFile(const label comm, const fileName& name);


    // Member Functions

    // Characteristics

        //- Communicator to use
        label comm() const noexcept
        {
            return comm_;
        }

        //- Set communicator to use [mutable]. Negative values are a no-op.
        //  \return old value
        label comm(label communicator) const noexcept
        {
            label old(comm_);
            if (communicator >= 0) comm_ = communicator;
            return old;
        }

        //- Distributed roots (parallel run)
        bool distributed() const noexcept
        {
            return distributed_;
        }

        //- Set distributed roots on/off [mutable]
        //  \return old value
        bool distributed(bool on) const noexcept
        {
            bool old(distributed_);
            distributed_ = on;
            return old;
        }

        //- The list of IO ranks (global ranks)
        //  Primarily for additional bookkeeping in non-parallel
        const labelList& ioRanks() const noexcept { return ioRanks_; }

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<fileOperation> info() const noexcept { return *this; }


    // Member Functions

        //- True if the fileOperation can be considered valid.
        //- At the moment, primarily used to detect the dummy fileOperation.
        virtual bool good() const { return true; }

        //- Transfer ownership of communicator to this fileOperation.
        //- Use with caution
        virtual void storeComm() const = 0;


        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const = 0;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const = 0;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return the file type: DIRECTORY, FILE or SYMLINK
            virtual fileName::Type type
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip = true,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip = true,
                const bool followLink = true
            ) const = 0;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true,
                const bool followLink = true
            ) const = 0;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const = 0;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const = 0;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const = 0;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const = 0;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const = 0;

            //- Remove a directory and its contents
            //  \param dir the directory to remove
            //  \param silent do not report missing directory
            //  \param emptyOnly only remove empty directories (recursive)
            virtual bool rmDir
            (
                const fileName& dir,
                const bool silent = false,
                const bool emptyOnly = false
            ) const = 0;

            //- Read dir/file (recursively if necessary) on master of the
            //- communicator, send and write contents to all 'writeOnProc'
            //- processors with local file name
            //  \param comm the communicator for broadcasting
            //  \param writeOnProc write on the processor
            //  \param src the source file/directory
            //  \param dst the target file/directory.
            //      If empty, treat as being identical to the src.
            virtual bool broadcastCopy
            (
                const label comm,
                const bool writeOnProc,
                const fileName& src,
                const fileName& dst
                // always recreates links
            ) const;


        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            //  Optional wanted typeName.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            //  Optional wanted typeName.
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName,
                const bool search = true
            ) const = 0;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject& io,
                const bool search = true
            ) const = 0;

            //- Search directory for objects. Used in IOobjectList.
            virtual fileNameList readObjects
            (
                const objectRegistry& db,
                const fileName& instance,
                const fileName& local,
                word& newInstance
            ) const;

            //- Read object header from supplied file
            virtual bool readHeader
            (
                IOobject&,
                const fileName&,
                const word& typeName
            ) const = 0;

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool readOnProc = true
            ) const = 0;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstreamOption::streamFormat format,
                const word& typeName
            ) const = 0;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state. Default action is to write to
            //  the objectPath using writeData. If !writeOnProc the
            //  file does not need to be written (this is used e.g. to
            //  suppress empty local lagrangian data)
            virtual bool writeObject
            (
                const regIOobject& io,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;


        // Filename (not IOobject) operations

            //- Search for a file or directory.
            //- Use IOobject version in preference
            virtual fileName filePath
            (
                const fileName&,
                const bool checkGzip = true,
                const bool followLink = true
            ) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream(const fileName&) const = 0;

            //- Generate an OSstream that writes a file
            virtual autoPtr<OSstream> NewOFstream
            (
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const = 0;

            //- Generate an OSstream that writes a file
            virtual autoPtr<OSstream> NewOFstream
            (
                IOstreamOption::atomicType atomic,
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const = 0;


        // File modification checking

            //- Add watching of a file. Returns handle
            virtual label addWatch(const fileName&) const;

            //- Remove watch on a file (using handle)
            virtual bool removeWatch(const label) const;

            //- Find index (or -1) of file in list of handles
            virtual label findWatch
            (
                const labelList& watchIndices,
                const fileName&
            ) const;

            //- Helper: add watches for list of regIOobjects
            virtual void addWatches(regIOobject&, const fileNameList&) const;

            //- Get name of file being watched (using handle)
            virtual fileName getFile(const label) const;

            //- Update state of all files
            virtual void updateStates
            (
                const bool masterOnly,
                const bool syncPar
            ) const;

            //- Get current state of file (using handle)
            virtual fileMonitor::fileState getState(const label) const;

            //- Set current state of file (using handle) to unmodified
            virtual void setUnmodified(const label) const;


        // Other

            //- Actual name of processors dir
            //- (for use in mode PROCOBJECT, PROCINSTANCE)
            virtual word processorsDir(const IOobject& io) const
            {
                return processorsBaseDir;
            }

            //- Actual name of processors dir
            //- (for use in mode PROCOBJECT, PROCINSTANCE)
            virtual word processorsDir(const fileName&) const
            {
                return processorsBaseDir;
            }

            //- Overall number of processors,
            //- from UPstream::nProcs() or detected from directories/results.
            label nProcs() const noexcept { return nProcs_; }

            //- Set number of processor directories/results.
            //  Used to cache format of e.g. processorsDDD.
            //  Returns old number of processors.
            //  Only used in decomposePar
            label nProcs(const label numProcs) noexcept
            {
                label old(nProcs_);
                nProcs_ = numProcs;
                return old;
            }

            //- Get number of processor directories/results.
            //  Used for e.g. reconstructPar, argList checking
            virtual label nProcs
            (
                const fileName& dir,
                const fileName& local = "",
                //! >0 if we know in advance what we're looking for
                const label wantedNProcs = 0
            ) const;

            //- Get sorted list of times
            virtual instantList findTimes(const fileName&, const word&) const;

            //- Find time instance where IOobject is located.
            //- The name of the IOobject can be empty, in which case only the
            //- IOobject::local() is checked.
            //- Does not search beyond \c stopInstance (if set) or \c constant.
            //  If the instance cannot be found:
            //  - FatalError when readOpt is (MUST_READ or READ_MODIFIED)
            //  - returns the \c stopInstance (if set and reached)
            //  - return \c constant if constant_fallback is true.
            //  - return an empty word if constant_fallback is false.
            //  .
            virtual IOobject findInstance
            (
                const IOobject& io,
                const scalar startValue,
                //! The search stop instance
                const word& stopInstance,
                //! Return \c "constant" instead of \c "" if the search failed
                const bool constant_fallback = true
            ) const;

            //- Callback for time change
            virtual void setTime(const Time&) const
            {}

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Forcibly parallel sync
            virtual void sync();

            //- Generate path (like io.path) from root+casename with any
            //  'processorXXX' replaced by procDir (usually 'processsors')
            fileName processorsCasePath
            (
                const IOobject& io,
                const word& procDir
            ) const;

            //- Generate path (like io.path) with provided instance and any
            //  'processorXXX' replaced by procDir (usually 'processsors')
            fileName processorsPath
            (
                const IOobject& io,
                const word& instance,
                const word& procDir
            ) const;

            //- Operating on fileName: replace processorXXX with procDir
            fileName processorsPath(const fileName&, const word& procDir) const;

            //- Split objectPath into part before 'processor' and part after.
            //
            //  Returns -1 or processor number and optionally number
            //  of processors. Use with care.
            //  - path/"processor"+Foam::name(proci)/local reconstructs input
            //  - path/"processors"+Foam::name(nProcs)/local reconstructs
            //    collated processors equivalence
            static label splitProcessorPath
            (
                const fileName& objectPath,
                fileName& path,
                fileName& procDir,
                fileName& local,
                procRangeType& group,
                label& nProcs
            );

            //- Detect processor number from '/aa/bb/processorDDD/cc'
            static label detectProcessorPath(const fileName& objPath);


    // Rank selection/sub-selection

        //- Get (contiguous) range/bounds of ranks addressed
        //- within the given main io-ranks.
        static labelRange subRanks(const labelUList& mainIOranks);

        //- Get list of global IO master ranks based on the hostname.
        //- It is assumed that each host range is contiguous.
        static labelList getGlobalHostIORanks();

        //- Get list of global IO ranks from FOAM_IORANKS env variable.
        //- If set, these correspond to the IO master ranks.
        static labelList getGlobalIORanks();


    // Housekeeping

        //- Same as nProcs
        label setNProcs(label numProcs) { return nProcs(numProcs); }
};


//- Read pathType as an integer value
inline Istream& operator>>(Istream& is, fileOperation::pathType& b)
{
    int val(0);
    is >> val;

    b = static_cast<fileOperation::pathType>(val);
    return is;
}

//- Write pathType as an integer value
inline Ostream& operator<<(Ostream& os, const fileOperation::pathType b)
{
    os << static_cast<int>(b);
    return os;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the current file handler
//- (will create default file handler if necessary).
//- Forwards to fileOperation::handler()
inline const fileOperation& fileHandler()
{
    return fileOperation::fileHandler();
}

//- Delete current file handler - forwards to fileOperation::handler()
//  Should have [[nodiscard]], but gcc ignores void casting.
inline refPtr<fileOperation> fileHandler(std::nullptr_t)
{
    return fileOperation::fileHandler(nullptr);
}

//- Replace the current file handler - forwards to fileOperation::handler().
//  \note legacy behaviour, so returns autoPtr instead of refPtr!
autoPtr<fileOperation> fileHandler(autoPtr<fileOperation>&& newHandler);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
