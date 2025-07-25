/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::fileOperations::uncollatedFileOperation

Description
    fileOperation that assumes file operations are local.

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_uncollatedFileOperation_H
#define Foam_fileOperations_uncollatedFileOperation_H

#include "fileOperation.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                     Class uncollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class uncollatedFileOperation
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

    // Protected Member Functions

        //- Search for an object.
        //    checkGlobal : also check undecomposed case
        //    isFile      : true:check for file  false:check for directory
        fileName filePathInfo
        (
            const bool checkGlobal,
            const bool isFile,
            const IOobject& io,
            const bool search
        ) const;

        //- Lookup name of processorsDDD using cache.
        //  \note Do not use any parallel synchronisation
        //  \return empty fileName if not found.
        virtual refPtr<dirIndexList> lookupProcessorsPath
        (
            const fileName& objectPath
        ) const;


public:

    //- Runtime type information
    TypeName("uncollated");


    // Constructors

        //- Default construct
        explicit uncollatedFileOperation(bool verbose = false);

        //- Construct from communicator with specified io-ranks
        explicit uncollatedFileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


    //- Destructor
    virtual ~uncollatedFileOperation();


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


        // (reg)IOobject functionality

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
            //  to read the contents.
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
                IOstreamOption::atomicType,
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
