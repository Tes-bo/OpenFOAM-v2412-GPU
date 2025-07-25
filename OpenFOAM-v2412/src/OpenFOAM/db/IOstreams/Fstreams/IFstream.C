/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "IFstream.H"
#include "OSspecific.H"  // For isFile(), fileSize()

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IFstream, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::DynamicList<char>
Foam::IFstream::readContents(IFstream& ifs)
{
    DynamicList<char> buffer;

    const auto inputSize = ifs.fileSize();

    if (inputSize <= 0)
    {
        // Nothing to read
    }
    else if (IOstreamOption::COMPRESSED == ifs.compression())
    {
        auto& iss = ifs.stdStream();

        // For compressed files, no idea how large the result will be.
        // So read chunk-wise.
        // Using the compressed size for the chunk size:
        // 50% compression = 2 iterations
        // 66% compression = 3 iterations
        // ...

        const uint64_t chunkSize =
        (
            (inputSize <= 1024)
          ? uint64_t(4096)
          : uint64_t(2*inputSize)
        );

        uint64_t beg = 0;

        for (int iter = 1; iter < 100000; ++iter)
        {
            // Manual resizing to use incremental vs doubling
            buffer.setCapacity(label(iter * chunkSize));
            buffer.resize(buffer.capacity());

            ifs.readRaw(buffer.data() + beg, chunkSize);
            const std::streamsize nread = iss.gcount();

            if
            (
                nread < 0
             || nread == std::numeric_limits<std::streamsize>::max()
            )
            {
                // Failed, but treat as normal 'done'
                buffer.resize(label(beg));
                break;
            }
            else
            {
                beg += uint64_t(nread);
                if (nread >= 0 && uint64_t(nread) < chunkSize)
                {
                    // normalExit = true;
                    buffer.resize(label(beg));
                    break;
                }
            }
        }
    }
    else
    {
        // UNCOMPRESSED
        {
            auto& iss = ifs.stdStream();

            buffer.setCapacity(label(inputSize));
            buffer.resize(buffer.capacity());

            ifs.readRaw(buffer.data(), buffer.size_bytes());
            const std::streamsize nread = iss.gcount();

            if
            (
                nread < 0
             || nread == std::numeric_limits<std::streamsize>::max()
            )
            {
                // Failed, but treat as normal 'done'
                buffer.clear();
            }
            else
            {
                buffer.resize(label(nread));  // Safety
            }
        }
    }

    return buffer;
}


Foam::DynamicList<char>
Foam::IFstream::readContents(const fileName& pathname)
{
    if (!pathname.empty())
    {
        IFstream ifs(pathname, IOstreamOption::BINARY);

        if (ifs.good())
        {
            return readContents(ifs);
        }
    }

    return DynamicList<char>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::IFstream
(
    const fileName& pathname,
    IOstreamOption streamOpt
)
:
    Foam::ifstreamPointer(pathname, streamOpt),
    ISstream(*(ifstreamPointer::get()), pathname, streamOpt)
{
    IOstreamOption::compression(ifstreamPointer::whichCompression());

    setClosed();

    setState(ifstreamPointer::get()->rdstate());

    if (good())
    {
        setOpened();
    }
    else
    {
        setBad();
    }

    lineNumber_ = 1;

    if (debug)
    {
        if (pathname.empty())
        {
            InfoInFunction
                << "Cannot open empty file name"
                << Foam::endl;
        }
        else if (IOstreamOption::COMPRESSED == IOstreamOption::compression())
        {
            InfoInFunction
                << "Decompressing " << (this->name() + ".gz") << Foam::endl;
        }

        if (!opened())
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for input\n" << info() << Foam::endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::streamsize Foam::IFstream::fileSize() const
{
    const std::istream* ptr = ifstreamPointer::get();

    if (!ptr || this->name().empty())
    {
        return std::streamsize(-1);
    }

    off_t fileLen = -1;

    if (IOstreamOption::COMPRESSED == ifstreamPointer::whichCompression())
    {
        fileLen = Foam::fileSize(this->name() + ".gz");
    }
    else
    {
        // TBD: special handing for wrapped icharstream
        // if
        // (
        //     const Foam::icharstream* charstr
        //   = dynamic_cast<const Foam::icharstream*>(ptr)>(ptr)
        // )
        // {
        //     return charstr->capacity();
        // }

        fileLen = Foam::fileSize(this->name());
    }

    if (fileLen >= 0)
    {
        return std::streamsize(fileLen);
    }

    return std::streamsize(-1);
}


std::istream& Foam::IFstream::stdStream()
{
    std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


const std::istream& Foam::IFstream::stdStream() const
{
    const std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


void Foam::IFstream::rewind()
{
    if (IOstreamOption::COMPRESSED == ifstreamPointer::whichCompression())
    {
        lineNumber_ = 1;  // Reset line number
        ifstreamPointer::reopen_gz(this->name());
        setState(ifstreamPointer::get()->rdstate());
    }
    else
    {
        ISstream::rewind();
    }
}


void Foam::IFstream::print(Ostream& os) const
{
    os  << "IFstream: ";
    ISstream::print(os);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::IFstream& Foam::IFstream::operator()() const
{
    if (!good())
    {
        // Also checks .gz file
        if (Foam::isFile(this->name(), true))
        {
            check(FUNCTION_NAME);
            FatalIOError.exit();
        }
        else
        {
            FatalIOErrorInFunction(*this)
                << "File " << this->name() << " does not exist"
                << exit(FatalIOError);
        }
    }

    return const_cast<IFstream&>(*this);
}


// ************************************************************************* //
