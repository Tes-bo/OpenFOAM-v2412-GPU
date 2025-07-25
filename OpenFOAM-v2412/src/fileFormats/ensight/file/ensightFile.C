/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ensightFile.H"
#include "ensightReadFile.H"
#include "error.H"
#include "List.H"
#include <cstring>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef_ = false;

float Foam::ensightFile::undefValue_ = Foam::floatScalarVGREAT;

const char* const Foam::ensightFile::coordinates = "coordinates";


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Put integers, floats etc in binary or ascii.
template<class Type>
static inline void putPrimitive
(
    const Type& value,
    OFstream& os,
    const int fieldWidth
)
{
    auto& oss = os.stdStream();

    if (os.format() == IOstreamOption::BINARY)
    {
        oss.write(reinterpret_cast<const char*>(&value), sizeof(Type));
    }
    else
    {
        oss.width(fieldWidth);
        oss << value;
    }
    os.syncState();
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::ensightFile::hasUndef(const UList<float>& field)
{
    for (const float val : field)
    {
        if (std::isnan(val))
        {
            return true;
        }
    }

    return true;
}


bool Foam::ensightFile::hasUndef(const UList<double>& field)
{
    for (const double val : field)
    {
        if (std::isnan(val))
        {
            return true;
        }
    }

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightFile::init()
{
    // The ASCII formatting specs for ensight files
    setf
    (
        std::ios_base::scientific,
        std::ios_base::floatfield
    );
    precision(5);

    // Handle transient single-file timestep information
    auto& oss = OFstream::stdStream();

    if (OFstream::is_appending())
    {
        // Already positioned at the EOF (in append mode), but be certain
        oss.seekp(0, std::ios_base::end);
        origFileSize_ = oss.tellp();
    }
    else
    {
        origFileSize_ = 0;
    }

    int64_t begin_footer(-1);
    List<int64_t> offsets;

    if (OFstream::is_appending())
    {
        // Temporarily open for reading as well.
        // No race condition since no writing is done concurrently with the
        // reading
        IFstream is(OFstream::name(), OFstream::format());

        begin_footer =
            ensightReadFile::getTimeStepFooter
            (
                is,
                offsets
            );
    }

    timeStepOffsets_ = std::move(offsets);

    if (OFstream::is_appending() && begin_footer > 0)
    {
        oss.seekp(begin_footer);
        OFstream::syncState();
    }

    // InfoErr << "output at: " << label(begin_footer) << nl;
    // InfoErr
    //     << "footer: " << label(begin_footer)
    //     << " time-steps: " << offsets.size() << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFile::ensightFile
(
    std::nullptr_t,  // dispatch tag
    IOstreamOption::appendType append,
    const fileName& pathname,
    IOstreamOption::streamFormat fmt
)
:
    OFstream
    (
        (
            // Only use atomic when not appending
            (append == IOstreamOption::NO_APPEND)
          ? IOstreamOption::ATOMIC
          : IOstreamOption::NON_ATOMIC
        ),
        pathname,
        fmt,
        (
            // Change APPEND_APP -> APPEND_ATE (file rewriting)
            (append == IOstreamOption::APPEND_APP)
          ? IOstreamOption::APPEND_ATE
          : append
        )
    ),
    origFileSize_(0)
{
    init();
}


Foam::ensightFile::ensightFile
(
    IOstreamOption::appendType append,
    const fileName& pathname,
    IOstreamOption::streamFormat fmt
)
:
    ensightFile
    (
        nullptr,
        append,
        ensight::FileName(pathname),
        fmt
    )
{}


Foam::ensightFile::ensightFile
(
    IOstreamOption::appendType append,
    const fileName& path,
    const fileName& name,
    IOstreamOption::streamFormat fmt
)
:
    ensightFile
    (
        nullptr,
        append,
        path/ensight::FileName(name),
        fmt
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightFile::~ensightFile()
{
    (void) writeTimeStepFooter();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef() noexcept
{
    return allowUndef_;
}


// float Foam::ensightFile::undefValue() noexcept
// {
//     return undefValue_;
// }


bool Foam::ensightFile::allowUndef(bool on) noexcept
{
    bool old = allowUndef_;
    allowUndef_ = on;
    return old;
}


float Foam::ensightFile::undefValue(float value) noexcept
{
    // enable its use too
    allowUndef_ = true;

    float old = undefValue_;
    undefValue_ = value;
    return old;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightFile::writeString(const char* str, size_t len)
{
    // Output 79 chars (ASCII) or 80 chars (BINARY)
    char buf[80];
    if (len > 80) len = 80;

    // TBD: truncate at newline? (shouldn't really occur anyhow)

    std::copy_n(str, len, buf);
    std::fill_n(buf + len, (80 - len), '\0');  // Pad trailing with nul

    auto& oss = stdStream();

    if (format() == IOstreamOption::BINARY)
    {
        oss.write(buf, 80);
    }
    else
    {
        buf[79] = 0;  // Max 79 in ASCII

        // TBD: Extra safety - trap newline in ASCII?
        // char* p = ::strchr(buf, '\n');
        // if (p) *p = 0;

        oss << buf;
    }

    syncState();
}


void Foam::ensightFile::writeString(const char* str)
{
    writeString(str, strlen(str));
}


void Foam::ensightFile::writeString(const std::string& str)
{
    writeString(str.data(), str.size());
}


Foam::Ostream& Foam::ensightFile::write(const char* str)
{
    writeString(str, strlen(str));
    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const word& str)
{
    writeString(str.data(), str.size());
    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const std::string& str)
{
    writeString(str.data(), str.size());
    return *this;
}


// Same as OFstream::writeRaw(buf, count)
Foam::Ostream& Foam::ensightFile::write
(
    const char* buf,
    std::streamsize count
)
{
    stdStream().write(buf, count);
    syncState();
    return *this;
}


void Foam::ensightFile::writeInt(const int32_t val, const int fieldWidth)
{
    putPrimitive<int32_t>(val, *this, fieldWidth);
}


void Foam::ensightFile::writeInt(const int64_t val, const int fieldWidth)
{
    putPrimitive<int32_t>(narrowInt32(val), *this, fieldWidth);
}


void Foam::ensightFile::writeFloat(const float val, const int fieldWidth)
{
    putPrimitive<float>(val, *this, fieldWidth);
}


void Foam::ensightFile::writeFloat(const double val, const int fieldWidth)
{
    putPrimitive<float>(narrowFloat(val), *this, fieldWidth);
}


Foam::Ostream& Foam::ensightFile::write(const int32_t val)
{
    putPrimitive<int32_t>(val, *this, 10);

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const int64_t val)
{
    putPrimitive<int32_t>(narrowInt32(val), *this, 10);

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const float val)
{
    putPrimitive<float>(val, *this, 12);

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const double val)
{
    putPrimitive<float>(narrowFloat(val), *this, 12);

    return *this;
}


void Foam::ensightFile::newline()
{
    if (format() == IOstreamOption::ASCII)
    {
        OFstream::write('\n');
    }
}


void Foam::ensightFile::writeUndef()
{
    write(undefValue_);
}


Foam::Ostream& Foam::ensightFile::writeKeyword(const keyType& key)
{
    if (allowUndef_)
    {
        writeString(key + " undef");
        newline();
        write(undefValue_);
        newline();
    }
    else
    {
        writeString(key);
        newline();
    }

    return *this;
}


void Foam::ensightFile::writeBinaryHeader()
{
    if (format() == IOstreamOption::BINARY)
    {
        writeString("C Binary");
        // newline();  // A no-op in binary
    }
}


// Footer information looks like this
//
/*  |---------------|---------------|-----------------------|
 *  | ASCII         | BINARY        | element               |
 *  |---------------|---------------|-----------------------|
 *  | "%20lld\n"    | int32         | nSteps                |
 *  | "%20lld\n"    | int64         | offset step 1         |
 *  | "%20lld\n"    | int64         | offset step 2         |
 *  | "%20lld\n"    | ..            |                       |
 *  | "%20lld\n"    | int64         | offset step n         |
 *  | "%20lld\n"    | int32         | flag (unused)         |
 *  | "%20lld\n"    | int64         | offset to nSteps      |
 *  | "%s\n"        | char[80]      | 'FILE_INDEX'          |
 *  |---------------|---------------|-----------------------|
 */

int64_t Foam::ensightFile::writeTimeStepFooter()
{
    if (timeStepOffsets_.empty())
    {
        return -1;
    }

    auto& oss = OFstream::stdStream();

    // The footer begin, which is also the current position
    const int64_t footer_begin(oss.tellp());

    // nSteps
    putPrimitive<int32_t>(int32_t(timeStepOffsets_.size()), *this, 20);
    newline();

    // offset step 1, 2, ... N
    for (int64_t off : timeStepOffsets_)
    {
        putPrimitive<int64_t>(off, *this, 20);
        newline();
    }

    // flag (unused)
    putPrimitive<int32_t>(0, *this, 20);
    newline();

    // The footer begin == position of nSteps
    putPrimitive<int64_t>(footer_begin, *this, 20);
    newline();

    // FILE_INDEX is "%s\n", not "%79s\n"
    // but our ASCII strings are truncated (nul-padded) anyhow

    writeString("FILE_INDEX");
    newline();

    // Reposition to begin of footer so that any subsequent output
    // will overwrite the footer too
    oss.seekp(footer_begin);

    return footer_begin;
}


//
// Convenience Output Methods
//

int64_t Foam::ensightFile::beginTimeStep()
{
    writeString("BEGIN TIME STEP");
    newline();

    auto& oss = OFstream::stdStream();

    const int64_t curr_pos(oss.tellp());
    timeStepOffsets_.push_back(curr_pos);

    // To avoid partly incomplete/incorrect footer information,
    // overwrite original footer if needed.

    if (curr_pos >= 0 && curr_pos < origFileSize_)
    {
        const char fill[] = "deadbeef";

        for
        (
            int64_t pos = curr_pos;
            pos < origFileSize_ && bool(oss);
            pos += 8
        )
        {
            // Overwrite with specified "junk" to avoid/detect corrupt
            // files etc. Don't worry about slightly increasing the
            // file size (ie, max 7 bytes) - it's unimportant
            oss.write(fill, 8);
        }

        // Maintain the original output position
        oss.seekp(curr_pos);

        OFstream::syncState();
    }

    return curr_pos;
}


int64_t Foam::ensightFile::endTimeStep()
{
    writeString("END TIME STEP");
    newline();

    return int64_t(stdStream().tellp());
}


void Foam::ensightFile::beginPart(const label index)
{
    writeString("part");
    newline();
    write(index+1); // Ensight starts with 1
    newline();
}


void Foam::ensightFile::beginPart
(
    const label index,
    const std::string& description
)
{
    beginPart(index);
    writeString(description);
    newline();
}


void Foam::ensightFile::beginCoordinates(const label npoints)
{
    writeString("coordinates");
    newline();

    write(npoints);
    newline();
}


void Foam::ensightFile::beginParticleCoordinates(const label nparticles)
{
    writeString("particle coordinates");
    newline();
    writeInt(nparticles, 8);  // Warning: unusual width
    newline();
}


void Foam::ensightFile::writeLabels(const UList<label>& list)
{
    for (const label val : list)
    {
        write(val);
        newline();
    }
}


void Foam::ensightFile::writeList(const UList<label>& field)
{
    for (const label val : field)
    {
        write(float(val));
        newline();
    }
}


void Foam::ensightFile::writeList(const UList<float>& field)
{
    for (const float val : field)
    {
        if (std::isnan(val))
        {
            writeUndef();
        }
        else
        {
            write(val);
        }
        newline();
    }
}


void Foam::ensightFile::writeList(const UList<double>& field)
{
    for (const double val : field)
    {
        if (std::isnan(val))
        {
            writeUndef();
        }
        else
        {
            write(val);
        }
        newline();
    }
}


// ************************************************************************* //
