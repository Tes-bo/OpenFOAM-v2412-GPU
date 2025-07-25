/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
#include "Pstream.H"
#include "OSspecific.H"
#include <fstream>
#include <cinttypes>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Implementation for broadcasting an individual file
static void broadcastFile_single
(
    const label comm,
    const bool writeOnProc,
    const fileName& srcName,
    const fileName& dstName,
    DynamicList<char>& buffer
)
{
    if (fileOperation::debug)
    {
        InfoInFunction
            << "Reading file " << srcName
            << " on master processor and copying to " << dstName
            << endl;
    }

    // Read file on master, broadcast to all but write on IOranks only.
    // This is a lot easier / possibly quicker?
    // than -allocating our own IO communicator
    // -send to all IO subProcs -releasing communicator

    // The file length uses uint64_t instead of off_t
    // since MINGW has something odd for off_t:
    //
    // ... error: call of overloaded 'min(off_t&, off_t&)' is ambiguous

    // Tuple2<off_t, mode_t> lengthAndMode(0, 0);
    Tuple2<uint64_t, mode_t> lengthAndMode(0, 0);

    std::unique_ptr<std::ifstream> srcStream;
    std::unique_ptr<std::ofstream> dstStream;

    if (UPstream::master(comm))
    {
        // Read (see newIFstream)
        auto fileLen = Foam::fileSize(srcName);
        if (fileLen > 0)
        {
            lengthAndMode.first() = uint64_t(fileLen);
        }
        lengthAndMode.second() = Foam::mode(srcName);

        srcStream.reset
        (
            new std::ifstream
            (
                srcName,
                std::ios_base::in | std::ios_base::binary
            )
        );
        if (!srcStream->good())
        {
            FatalIOErrorInFunction(srcName)
                << "Could not open file for reading!"
                << exit(FatalIOError);
        }
    }

    if
    (
        writeOnProc
     &&
        (
            UPstream::master(comm)
          ? (srcName != dstName)
          : UPstream::is_subrank(comm)
        )
    )
    {
        // Make sure the destination directory exists.
        // - will fail itself if not possible
        // - no-op if directory already exists
        Foam::mkDir(dstName.path());

        dstStream.reset
        (
            new std::ofstream
            (
                dstName,
                std::ios_base::out | std::ios_base::binary
            )
        );

        if (!dstStream->good())
        {
            // Fail noisily or silently?
            if (!dstStream->good())
            {
                dstStream.reset(nullptr);
            }
        }

        // Adjust mode?
        // Foam::chMode(dstName, fileMode);
    }

    // Broadcast size and mode (contiguous data)
    UPstream::broadcast
    (
        reinterpret_cast<char*>(&lengthAndMode),
        sizeof(lengthAndMode),
        comm
    );

    uint64_t fileLength = lengthAndMode.first();

    const uint64_t maxChunkSize =
    (
        (UPstream::maxCommsSize > 0)
      ? uint64_t(UPstream::maxCommsSize)
      : (UPstream::maxCommsSize < 0)    // (numBytes fewer than INT_MAX)
      ? uint64_t(INT_MAX + UPstream::maxCommsSize)
      : uint64_t(INT_MAX)               // MPI limit is <int>
    );


    while (fileLength > 0)
    {
        const uint64_t sendSize = std::min(fileLength, maxChunkSize);
        fileLength -= sendSize;

        // Read file contents into a character buffer
        buffer.resize_nocopy(static_cast<label>(sendSize));

        if (srcStream)
        {
            srcStream->read(buffer.data_bytes(), buffer.size_bytes());
        }

        UPstream::broadcast(buffer.data_bytes(), buffer.size_bytes(), comm);

        if (dstStream)
        {
            dstStream->write(buffer.data_bytes(), buffer.size_bytes());
        }
    }
}


// Implementation for broadcasting directory contents or an individual file
static bool broadcastFile_recursive
(
    const label comm,
    const bool writeOnProc,
    const fileName& src,
    const fileName& dst,
    // const bool followLink
    DynamicList<char>& buffer
)
{
    // Read file on master, broadcast to all but write on procs with doWrite
    // only. This is a lot easier / possibly quicker?
    // than -allocating our own IO communicator
    // -send to all IO subProcs -releasing communicator

    if (fileOperation::debug)
    {
        InfoInFunction
            << "Reading " << src
            << " on master processor and writing a copy to " << dst
            << endl;
    }

    fileName::Type srcType;
    if (UPstream::master(comm))
    {
        srcType = src.type(false);  // followLink = false, gzip = false
    }

    UPstream::broadcast
    (
        reinterpret_cast<char*>(&srcType),  // contiguous data
        sizeof(srcType),
        comm
    );

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        broadcastFile_single(comm, writeOnProc, src, dst, buffer);
    }
    else if (srcType == fileName::SYMLINK)
    {
        WarningInFunction<< "Copying symbolic links not supported" << endl;

        // Read the link target
        fileName linkTarget;
        if (UPstream::master(comm))
        {
            linkTarget = Foam::readLink(src);
        }
        Pstream::broadcast(linkTarget, comm);

        if
        (
            writeOnProc
         &&
            (
                UPstream::master(comm)
              ? (src != dst)
              : UPstream::is_subrank(comm)
            )
        )
        {
            // Recreate softlink on remote processor
            return Foam::ln(linkTarget, dst);
        }
    }
    else if (srcType == fileName::DIRECTORY)
    {
        // Copy files
        {
            fileNameList files;
            if (UPstream::master(comm))
            {
                files = Foam::readDir
                (
                    src,
                    fileName::FILE,
                    false  // Never trim '.gz' ending
                    //followLink
                );
            }
            Pstream::broadcast(files, comm);

            for (const fileName& item : files)
            {
                // File to file
                broadcastFile_recursive
                (
                    comm,
                    writeOnProc,
                    src/item,
                    dst/item,
                    //followLink
                    buffer
                );
            }
        }

        // Copy softlinks
        {
            fileNameList files;
            if (UPstream::master(comm))
            {
                files = Foam::readDir
                (
                    src,
                    fileName::SYMLINK,
                    false  // Never trim '.gz' ending
                    //followLink
                );
            }
            Pstream::broadcast(files, comm);

            for (const fileName& item : files)
            {
                // Softlink to softlink
                broadcastFile_recursive
                (
                    comm,
                    writeOnProc,
                    src/item,
                    dst/item,
                    //followLink
                    buffer
                );
            }
        }


        // Copy sub directories
        {
            fileNameList dirs;
            if (UPstream::master(comm))
            {
                dirs = Foam::readDir
                (
                    src,
                    fileName::DIRECTORY,
                    false  // Never trim '.gz' ending
                    //followLink
                );
            }
            Pstream::broadcast(dirs, comm);

            for (const fileName& item : dirs)
            {
                // Dir to Dir
                broadcastFile_recursive
                (
                    comm,
                    writeOnProc,
                    src/item,
                    dst/item,
                    //followLink
                    buffer
                );
            }
        }
    }
    else if (srcType == fileName::UNDEFINED)
    {
        WarningInFunction<< "No known file type: " << src << endl;
        return false;
    }
    else
    {
        return false;
    }

    return true;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::fileOperation::broadcastCopy
(
    const label comm,
    const bool writeOnProc,
    const fileName& srcPath,
    const fileName& dstPath
    //const bool followLink
) const
{
    DynamicList<char> fileContents;


    // FUTURE:
    // Handling with broadcast of file -> directory etc
    // as per Foam::cp()
    //
    // However, this adds extra communicattion and logic etc...

#if 0
    const fileName::Type srcType = src.type(followLink);

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }
    }
    else if (srcType == fileName::DIRECTORY)
    {
        if (destFile.type() == fileName::DIRECTORY)
        {
            // Both are directories. Could mean copy contents or copy
            // recursively.  Don't actually know what the user wants,
            // but assume that if names are identical == copy contents.
            //
            // So: "path1/foo" "path2/foo"  copy contents
            // So: "path1/foo" "path2/bar"  copy directory

            const word srcDirName = src.name();
            if (destFile.name() != srcDirName)
            {
                destFile /= srcDirName;
            }
        }
    }
#endif

    return broadcastFile_recursive
    (
        comm,
        writeOnProc,
        srcPath,
        (dstPath.empty() ? srcPath : dstPath),
        fileContents
    );
}


// ************************************************************************* //
