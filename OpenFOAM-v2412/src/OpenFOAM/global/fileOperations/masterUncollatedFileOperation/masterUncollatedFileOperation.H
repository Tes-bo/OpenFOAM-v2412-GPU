/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

Class
    Foam::fileOperations::masterUncollatedFileOperation

Description
    fileOperations that performs all file operations on the master processor.
    Requires the calls to be parallel synchronised!

    Limitations
    - no /processor in filename
    - no /uniform/ in the filename

    The main logic is in ::filePath which returns a
    - same path on all processors. This can either be a global file
      (system/controlDict, processorXXX/0/uniform/) or a collated file
      (processors/0/p)
    - same path on all processors of the local communicator
      (processors4_0-1/0/p)
    - different path on all processors (processor0/0/p)

    system/controlDict:
    filePath worldmaster: \<globalRoot\>/system/controlDict
             localmaster:           ,,
             slave      :           ,,

    processor0/uniform/time
    filePath worldmaster: \<globalRoot\>/processorXXX/uniform/time
             localmaster:           ,,
             slave      :           ,,

    processors0/0/p
    processors10/0/p
    processors10_2-4/0/p

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_masterUncollatedFileOperation_H
#define Foam_fileOperations_masterUncollatedFileOperation_H

#include "fileOperation.H"
#include "OSspecific.H"
#include "HashPtrTable.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class PstreamBuffers;

namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                     Class masterUncollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class masterUncollatedFileOperation
:
    public fileOperation
{
    // Private Data

        //- Communicator allocated/managed by us
        mutable label managedComm_;


    // Private Member Functions

        //- Any initialisation steps after constructing
        void init(bool verbose);


protected:

    // Protected Data

        //- Cached times for a given directory
        mutable HashPtrTable<DynamicList<instant>> times_;


    // Protected Operation Functors

        class mkDirOp
        {
            const mode_t mode_;
        public:
            mkDirOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::mkDir(f, mode_);
            }
        };

        class chModOp
        {
            const mode_t mode_;
        public:
            chModOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::chMod(f, mode_);
            }
        };

        class modeOp
        {
            const bool followLink_;
        public:
            modeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            mode_t operator()(const fileName& f) const
            {
                return Foam::mode(f, followLink_);
            }
        };

        class typeOp
        {
            const bool followLink_;
        public:
            typeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            // Returns int (for reductions) instead of fileName::Type
            int operator()(const fileName& f) const
            {
                return Foam::type(f, followLink_);
            }
        };

        class existsOp
        {
            const bool checkGzip_;
            const bool followLink_;
        public:
            existsOp(const bool checkGzip, const bool followLink)
            :
                checkGzip_(checkGzip),
                followLink_(followLink)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::exists(f, checkGzip_, followLink_);
            }
        };

        class isDirOp
        {
            const bool followLink_;
        public:
            isDirOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::isDir(f, followLink_);
            }
        };

        class isFileOp
        {
            const bool checkGzip_;
            const bool followLink_;
        public:
            isFileOp(const bool checkGzip, const bool followLink)
            :
                checkGzip_(checkGzip),
                followLink_(followLink)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::isFile(f, checkGzip_, followLink_);
            }
        };

        class fileSizeOp
        {
            const bool followLink_;
        public:
            fileSizeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            off_t operator()(const fileName& f) const
            {
                return Foam::fileSize(f, followLink_);
            }
        };

        class lastModifiedOp
        {
            const bool followLink_;
        public:
            lastModifiedOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            time_t operator()(const fileName& f) const
            {
                return Foam::lastModified(f, followLink_);
            }
        };

        class highResLastModifiedOp
        {
            const bool followLink_;
        public:
            highResLastModifiedOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            double operator()(const fileName& f) const
            {
                return Foam::highResLastModified(f, followLink_);
            }
        };

        class mvBakOp
        {
            std::string ext_;
        public:
            mvBakOp(const std::string& ext)
            :
                ext_(ext)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::mvBak(f, ext_);
            }
        };

        class rmOp
        {
        public:
            bool operator()(const fileName& f) const
            {
                return Foam::rm(f);
            }
        };

        class rmDirOp
        {
            bool silent_;
            bool emptyOnly_;
        public:
            rmDirOp(bool silent=false, bool emptyOnly=false)
            :
                silent_(silent),
                emptyOnly_(emptyOnly)
            {}

            bool operator()(const fileName& f) const
            {
                return Foam::rmDir(f, silent_, emptyOnly_);
            }
        };

        class cpOp
        {
            const bool followLink_;
        public:
            cpOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::cp(src, dest, followLink_);
            }
        };

        class lnOp
        {
        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::ln(src, dest);
            }
        };

        class mvOp
        {
            const bool followLink_;
        public:
            mvOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::mv(src, dest, followLink_);
            }
        };

        class fileOrNullOp
        {
            const bool isFile_;
        public:
            fileOrNullOp(const bool isFile)
            :
                isFile_(isFile)
            {}

            const fileName& operator()(const fileName& f) const
            {
                return
                (
                    (isFile_ ? Foam::isFile(f) : Foam::isDir(f))
                  ? f
                  : fileName::null  // const reference, not temporary
                );
            }
        };

        class readDirOp
        {
            const fileName::Type type_;
            const bool filtergz_;
            const bool followLink_;
        public:
            readDirOp
            (
                const fileName::Type type,
                const bool filtergz,
                const bool followLink
            )
            :
                type_(type),
                filtergz_(filtergz),
                followLink_(followLink)
            {}

            fileNameList operator()(const fileName& f) const
            {
                return Foam::readDir(f, type_, filtergz_, followLink_);
            }
        };


    // Private Member Functions

        template<class Type, class FileOp>
        Type masterOp
        (
            const fileName& fName,
            const FileOp& fop,
            const int tag,
            const label comm
        ) const;

        template<class Type, class FileOp>
        Type masterOp
        (
            const fileName& src,
            const fileName& dest,
            const FileOp& fop,
            const int tag,
            const label comm
        ) const;

        //- Search (locally!) for object; return info on how it was found.
        //  Does not do any parallel communication.
        //    checkGlobal : also check undecomposed case
        //    isFile      : true:check for file  false:check for directory
        // searchType    : how was found
        // processorsDir : name of processor directory
        // instance      : instance
        virtual fileName filePathInfo
        (
            const bool checkGlobal,
            const bool isFile,
            const IOobject& io,
            const dirIndexList& pDirs,
            const bool search,
            pathType& searchType,
            word& processorsDir,
            word& instance
        ) const;

        //- Construct filePath
        fileName localObjectPath
        (
            const IOobject&,
            const pathType& searchType,
            const word& processorsDir,
            const word& instancePath
        ) const;

        //- Read file contents and send to processors.
        //  Handles compressed or uncompressed files
        static void readAndSend
        (
            const fileName& filePath,
            const labelUList& recvProcs,
            PstreamBuffers& pBufs
        );

        //- Read files on comms master
        static autoPtr<ISstream> read
        (
            IOobject& io,
            const label comm,
            const bool uniform,             // on comms master only
            const fileNameList& filePaths,  // on comms master and sub-ranks
            const boolUList& readOnProcs    // on comms master and sub-ranks
        );

        //- Helper: check IO for local existence. Like filePathInfo but
        //  without parent searchign and instance searching
        bool exists(const dirIndexList&, IOobject& io) const;


public:

        //- Runtime type information
        TypeName("masterUncollated");


    // Static Data

        //- Max size of parallel communications. Switches from non-blocking
        //  to scheduled when reading/writing files. Read as float to enable
        //  easy specification of large sizes.
        static float maxMasterFileBufferSize;


    // Constructors

        //- Default construct
        explicit masterUncollatedFileOperation(bool verbose = false);

        //- Construct from communicator with specified io-ranks
        explicit masterUncollatedFileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


    //- Destructor
    virtual ~masterUncollatedFileOperation();


    // Member Functions

        //- Transfer ownership of communicator to this fileOperation.
        //- Use with caution
        virtual void storeComm() const;


        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Return the file type: DIRECTORY, FILE or SYMLINK
            virtual fileName::Type type
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true,
                const bool followLink = true
            ) const;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const;

            //- Remove a directory and its contents
            //  \param dir the directory to remove
            //  \param silent do not report missing directory
            //  \param emptyOnly only remove empty directories (recursive)
            virtual bool rmDir
            (
                const fileName& dir,
                const bool silent = false,
                const bool emptyOnly = false
            ) const;


        // (reg)IOobject functinality

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject& io,
                const word& typeName,
                const bool search
            ) const;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject& io,
                const bool search
            ) const;

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
            ) const;

            //- Reads header for regIOobject and returns an ISstream
            //- to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool readOnProc = true
            ) const;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstreamOption::streamFormat format,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject& io,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream(const fileName&) const;

            //- Generate an OSstream that writes a file
            virtual autoPtr<OSstream> NewOFstream
            (
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;

            //- Generate an OSstream that writes a file
            virtual autoPtr<OSstream> NewOFstream
            (
                IOstreamOption::atomicType atomic,
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;


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
            virtual void setTime(const Time&) const;

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Forcibly parallel sync
            virtual void sync();

            //- Return cached times
            const HashPtrTable<DynamicList<instant>>& times() const noexcept
            {
                return times_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "masterUncollatedFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
