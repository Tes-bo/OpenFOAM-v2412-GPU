/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "fstreamPointer.H"
#include "OCountStream.H"
#include "OSspecific.H"
#include <cstdio>

// HAVE_LIBZ defined externally
// #define HAVE_LIBZ

#ifdef HAVE_LIBZ
#include "gzstream.h"
#endif /* HAVE_LIBZ */

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::ifstreamPointer::supports_gz() noexcept
{
    #ifdef HAVE_LIBZ
    return true;
    #else
    return false;
    #endif
}


bool Foam::ofstreamPointer::supports_gz() noexcept
{
    #ifdef HAVE_LIBZ
    return true;
    #else
    return false;
    #endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ifstreamPointer::ifstreamPointer
(
    const fileName& pathname,
    IOstreamOption streamOpt  // Currently unused
)
:
    ptr_()
{
    open(pathname, streamOpt);
}


Foam::ifstreamPointer::ifstreamPointer
(
    const fileName& pathname
)
:
    ptr_()
{
    open(pathname);
}


Foam::ofstreamPointer::ofstreamPointer() noexcept
:
    ptr_(),
    mode_(modeType::NONE)
{}


Foam::ofstreamPointer::ofstreamPointer(std::nullptr_t)
:
    ptr_(new Foam::ocountstream),
    mode_(modeType::NONE)
{}


Foam::ofstreamPointer::ofstreamPointer
(
    const fileName& pathname,
    IOstreamOption streamOpt,
    IOstreamOption::appendType append,
    bool atomic
)
:
    ptr_(),
    mode_(modeType::NONE)
{
    // Leave std::ios_base::trunc implicitly handled to make things
    // easier for append mode.

    std::ios_base::openmode openmode
    (
        std::ios_base::out | std::ios_base::binary
    );

    if (append == IOstreamOption::APPEND_APP)
    {
        openmode |= std::ios_base::app;

        // Cannot append to gzstream
        streamOpt.compression(IOstreamOption::UNCOMPRESSED);

        // Cannot use append + atomic operation, without lots of extra work
        atomic = false;
    }
    else if (append == IOstreamOption::APPEND_ATE)
    {
        // Handle an "append-like" mode by opening "r+b" and NOT as "ab"
        // - file already exists: Sets read position to start
        // - file does not exist: Error

        openmode |= std::ios_base::in;  // [SIC] - use read bit, not append!

        // Cannot append to gzstream
        streamOpt.compression(IOstreamOption::UNCOMPRESSED);

        // Cannot use append + atomic operation, without lots of extra work
        atomic = false;
    }


    // When opening new files, remove file variants out of the way.
    // Eg, opening "file1"
    // - remove old "file1.gz" (compressed)
    // - also remove old "file1" if it is a symlink and we are not appending
    //
    // Not writing into symlinked files avoids problems with symlinked
    // initial fields (eg, 0/U -> ../0.orig/U)

    const fileName pathname_gz(pathname + ".gz");
    const fileName pathname_tmp(pathname + "~tmp~");

    fileName::Type fType = fileName::Type::UNDEFINED;

    if (IOstreamOption::COMPRESSED == streamOpt.compression())
    {
        // Output compression requested.

        #ifdef HAVE_LIBZ
        // TBD:
        // atomic = true;  // Always treat COMPRESSED like an atomic

        const fileName& target = (atomic ? pathname_tmp : pathname_gz);

        // Remove old uncompressed version (if any)
        fType = Foam::type(pathname, false);
        if (fType == fileName::SYMLINK || fType == fileName::FILE)
        {
            Foam::rm(pathname);
        }

        // Avoid writing into symlinked files (non-append mode)
        if (atomic || (append == IOstreamOption::NO_APPEND))
        {
            fType = Foam::type(target, false);
            if (fType == fileName::SYMLINK)
            {
                Foam::rm(target);
            }
        }

        ptr_.reset(new ogzstream(target, openmode));

        #else /* HAVE_LIBZ */

        streamOpt.compression(IOstreamOption::UNCOMPRESSED);

        Warning
            << nl
            << "No write support for gz compressed files (libz)"
            << " : downgraded to UNCOMPRESSED" << nl
            << "file: " << pathname_gz << endl;

        #endif /* HAVE_LIBZ */
    }

    if (IOstreamOption::COMPRESSED != streamOpt.compression())
    {
        const fileName& target = (atomic ? pathname_tmp : pathname);

        // Remove old compressed version (if any)
        fType = Foam::type(pathname_gz, false);
        if (fType == fileName::SYMLINK || fType == fileName::FILE)
        {
            Foam::rm(pathname_gz);
        }

        // Avoid writing into symlinked files (non-append mode)
        if (atomic || (append == IOstreamOption::NO_APPEND))
        {
            fType = Foam::type(target, false);
            if (fType == fileName::SYMLINK)
            {
                Foam::rm(target);
            }
        }

        // File pointer (std::ofstream)
        auto filePtr = std::make_unique<std::ofstream>(target, openmode);

        if (append == IOstreamOption::APPEND_APP)
        {
            // Final handling for append 'app' (always non-atomic)

            // Set output position to the end (like std::ios_base::ate)
            // but only to test if the file had a size.
            // No real performance problem since any subsequent write
            // will do the same anyhow.

            filePtr->seekp(0, std::ios_base::end);
            if (filePtr->tellp() <= 0)
            {
                // Did not open an existing file
                append = IOstreamOption::NO_APPEND;
            }
        }
        else if (append == IOstreamOption::APPEND_ATE)
        {
            // Final handling for append 'ate' (always non-atomic)

            if (filePtr->good())
            {
                // Success if file already exists.
                // Set output position to the end - like std::ios_base::ate

                filePtr->seekp(0, std::ios_base::end);

                if (filePtr->tellp() <= 0)
                {
                    // Did not open an existing file
                    append = IOstreamOption::NO_APPEND;
                }
            }
            else
            {
                // Error if file does not already exist.
                // Reopen as regular output mode only.

                // Did not open an existing file
                append = IOstreamOption::NO_APPEND;

                if (filePtr->is_open())
                {
                    filePtr->close();
                }
                filePtr->clear();
                filePtr->open
                (
                    target,
                    (std::ios_base::out | std::ios_base::binary)
                );
            }
        }
        ptr_.reset(filePtr.release());
    }

    // Is appending to an existing file
    if (append != IOstreamOption::NO_APPEND)
    {
        mode_ = modeType::APPENDING;
    }

    // An atomic output operation (normally not appending!)
    if (atomic)
    {
        mode_ |= modeType::ATOMIC;
    }
}


Foam::ofstreamPointer::ofstreamPointer
(
    const fileName& pathname,
    IOstreamOption::compressionType comp,
    IOstreamOption::appendType append,
    bool atomic
)
:
    ofstreamPointer
    (
        pathname,
        IOstreamOption(IOstreamOption::ASCII, comp),
        append,
        atomic
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ifstreamPointer::open
(
    const fileName& pathname,
    IOstreamOption streamOpt  // Currently unused
)
{
    // Forcibly close old stream (if any)
    ptr_.reset(nullptr);

    const std::ios_base::openmode openmode
    (
        std::ios_base::in | std::ios_base::binary
    );

    ptr_.reset(new std::ifstream(pathname, openmode));

    if (!ptr_->good())
    {
        // Try compressed version instead

        const fileName pathname_gz(pathname + ".gz");

        if (Foam::isFile(pathname_gz, false))
        {
            #ifdef HAVE_LIBZ

            ptr_.reset(new igzstream(pathname_gz, openmode));

            #else /* HAVE_LIBZ */

            FatalError
                << "No read support for gz compressed files (libz)"
                << " : could use 'gunzip' from the command-line" << nl
                << "file: " << pathname_gz << endl
                << exit(FatalError);

            #endif /* HAVE_LIBZ */
        }
        else
        {
            // TBD:
            // Can also fallback and open .orig files too
            //
            // auto* file = dynamic_cast<std::ifstream*>(ptr_.get());
            // file->open(pathname + ".orig", mode);
        }
    }
}


void Foam::ifstreamPointer::reopen_gz(const std::string& pathname)
{
    #ifdef HAVE_LIBZ
    auto* gz = dynamic_cast<igzstream*>(ptr_.get());

    if (gz)
    {
        // Special treatment for gzstream
        gz->close();
        gz->clear();

        gz->open
        (
            pathname + ".gz",
            (std::ios_base::in | std::ios_base::binary)
        );
        return;
    }
    #endif /* HAVE_LIBZ */
}


void Foam::ofstreamPointer::reopen(const std::string& pathname)
{
    #ifdef HAVE_LIBZ
    auto* gz = dynamic_cast<ogzstream*>(ptr_.get());

    if (gz)
    {
        // Special treatment for gzstream
        gz->close();
        gz->clear();

        if (mode_ & modeType::ATOMIC)
        {
            gz->open
            (
                pathname + "~tmp~",
                (std::ios_base::out | std::ios_base::binary)
            );
        }
        else
        {
            gz->open
            (
                pathname + ".gz",
                (std::ios_base::out | std::ios_base::binary)
            );
        }
        return;
    }
    #endif /* HAVE_LIBZ */

    auto* file = dynamic_cast<std::ofstream*>(ptr_.get());

    if (file)
    {
        if (file->is_open())
        {
            file->close();
        }
        file->clear();

        // Invalidate the appending into existing file information
        // since rewind usually means overwrite
        mode_ &= ~modeType::APPENDING;

        if (mode_ & modeType::ATOMIC)
        {
            file->open
            (
                pathname + "~tmp~",
                (std::ios_base::out | std::ios_base::binary)
            );
        }
        else
        {
            file->open
            (
                pathname,
                (std::ios_base::out | std::ios_base::binary)
            );
        }
        return;
    }
}


void Foam::ofstreamPointer::close(const std::string& pathname)
{
    // Invalidate the appending into existing file information
    mode_ &= ~modeType::APPENDING;

    if (pathname.empty() || !(mode_ & modeType::ATOMIC))
    {
        return;
    }

    #ifdef HAVE_LIBZ
    auto* gz = dynamic_cast<ogzstream*>(ptr_.get());

    if (gz)
    {
        // Special treatment for gzstream
        gz->close();
        gz->clear();

        std::rename
        (
            (pathname + "~tmp~").c_str(),
            (pathname + ".gz").c_str()
        );
        return;
    }
    #endif /* HAVE_LIBZ */

    auto* file = dynamic_cast<std::ofstream*>(ptr_.get());

    if (file)
    {
        if (file->is_open())
        {
            file->close();
        }
        file->clear();

        std::rename
        (
            (pathname + "~tmp~").c_str(),
            pathname.c_str()
        );
        return;
    }
}


Foam::IOstreamOption::compressionType
Foam::ifstreamPointer::whichCompression() const
{
    #ifdef HAVE_LIBZ
    if (dynamic_cast<const igzstream*>(ptr_.get()))
    {
        return IOstreamOption::compressionType::COMPRESSED;
    }
    #endif /* HAVE_LIBZ */

    return IOstreamOption::compressionType::UNCOMPRESSED;
}


Foam::IOstreamOption::compressionType
Foam::ofstreamPointer::whichCompression() const
{
    #ifdef HAVE_LIBZ
    if (dynamic_cast<const ogzstream*>(ptr_.get()))
    {
        return IOstreamOption::compressionType::COMPRESSED;
    }
    #endif /* HAVE_LIBZ */

    return IOstreamOption::compressionType::UNCOMPRESSED;
}


// ************************************************************************* //
