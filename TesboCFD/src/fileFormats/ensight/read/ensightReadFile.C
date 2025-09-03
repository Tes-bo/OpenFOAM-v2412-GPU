/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ensightReadFile.H"
#include "stringOps.H"
#include "defineDebugSwitch.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineDebugSwitchWithName(Foam::ensightReadFile, "ensightReadFile", 0);

registerDebugSwitchWithName(Foam::ensightReadFile, ensight, "ensightReadFile");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Get integers, floats etc in binary or ascii.
template<class Type>
static inline Type getPrimitive(IFstream& is)
{
    Type value(0);

    auto& iss = is.stdStream();

    if (is.format() == IOstreamOption::BINARY)
    {
        iss.read(reinterpret_cast<char*>(&value), sizeof(Type));
    }
    else
    {
        iss >> value;
    }
    is.syncState();

    return value;
}


// Get an Ensight string value (binary or ascii).
static inline void readEnsightString(IFstream& is, std::string& value)
{
    if (is.format() == IOstreamOption::BINARY)
    {
        auto& iss = is.stdStream();

        // Binary string is *exactly* 80 characters
        value.resize(80, '\0');
        iss.read(&value[0], 80);
        const std::streamsize gcount = iss.gcount();
        value.erase(gcount <= 0 ? 0 : gcount);  // Truncated?

        // Could exit on truncated input, but no real advantage

        // Truncate at the first embedded '\0'
        const auto endp = value.find('\0');

        if (endp != std::string::npos)
        {
            value.erase(endp);
        }

        // May have been padded with trailing spaces - remove those
        stringOps::inplaceTrimRight(value);

        is.syncState();
    }
    else
    {
        value.clear();
        while (value.empty() && !is.eof())
        {
            is.getLine(value);
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

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

int64_t Foam::ensightReadFile::getTimeStepFooter
(
    IFstream& is,
    // File offsets for each time step (if any)
    List<int64_t>& offsets
)
{
    std::string buffer;

    auto& iss = is.stdStream();
    const auto lineNum = is.lineNumber();
    const auto curr_pos = iss.tellg();

    if (curr_pos < 0)
    {
        // Impossible positioning - exit
        is.lineNumber(lineNum);     // Restore line number
        offsets.clear();
        return -1;
    }

    iss.seekg(0, std::ios_base::end);
    const auto end_pos = iss.tellg();

    // As a minimum, expect at least 1 time step, so have four integers
    // (nSteps, offset step 1, flag, file offset) and the string (10 chars).
    // Thus always at least 80+ chars.

    if (end_pos <= 80)
    {
        // Looks quite impossible - exit

        is.lineNumber(lineNum);     // Restore line number
        iss.seekg(curr_pos);        // Restore file position

        offsets.clear();
        return -1;
    }

    // Get the last 80 chars as a character string
    iss.seekg(-80, std::ios_base::end);

    const auto fmt = is.format(IOstreamOption::BINARY);
    readEnsightString(is, buffer);
    is.format(fmt);

    int64_t footer_begin(0);

    const auto endp = buffer.find("FILE_INDEX");

    if (endp == std::string::npos)
    {
        // Not found

        is.lineNumber(lineNum);     // Restore line number
        iss.seekg(curr_pos);        // Restore file position

        offsets.clear();
        return -1;
    }
    else if (fmt == IOstreamOption::ASCII)
    {
        // In ASCII, the last 80 chars will also include a few integers
        buffer.erase(endp);  // Remove FILE_INDEX ...
        auto split = stringOps::splitSpace(buffer);

        if (!split.empty())
        {
            footer_begin = Foam::readInt64(split.back().str());
        }
    }
    else
    {
        // Position before string (80 bytes) and int64 value (8 bytes)
        iss.seekg(-88, std::ios_base::end);
        footer_begin = getPrimitive<int64_t>(is);
    }


    // The number of steps is stored as int32 at the beginning of the footer
    int32_t nSteps(0);

    if (footer_begin)
    {
        iss.seekg(footer_begin);
        nSteps = getPrimitive<int32_t>(is);
    }

    offsets.resize_nocopy(nSteps);

    // Next footer entries are the offsets per time-step
    for (int32_t step = 0; step < nSteps; ++step)
    {
        offsets[step] = getPrimitive<int64_t>(is);
    }

    is.lineNumber(lineNum);     // Restore line number
    iss.seekg(curr_pos);        // Restore file position

    return footer_begin;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightReadFile::readString(std::string& value)
{
    readEnsightString(*this, value);
}


void Foam::ensightReadFile::init(bool detectFormat)
{
    if (!IFstream::good())
    {
        FatalErrorInFunction
            << "Cannot read file " << IFstream::name() << nl
            << exit(FatalError);
    }

    auto& iss = stdStream();

    auto lineNum = lineNumber();
    auto curr_pos = iss.tellg();  // The starting position (should be 0)

    string buffer;

    if (detectFormat)
    {
        // Read initial string as BINARY
        IFstream::format(IOstreamOption::BINARY);

        readEnsightString(*this, buffer);

        // Detect BINARY vs ASCII by testing for initial "(C|Fortran) Binary"
        if (buffer.contains("Binary") || buffer.contains("binary"))
        {
            // Format is BINARY
            IFstream::format(IOstreamOption::BINARY);

            // New backtracking point is after the initial "C Binary" string
            curr_pos = iss.tellg();

            // Get the next (optional) line after the "C Binary" (if any)
            // and before the description.
            readEnsightString(*this, buffer);
        }
        else
        {
            // Not binary => ASCII
            IFstream::format(IOstreamOption::ASCII);

            // Rewind to the beginning again
            iss.seekg(curr_pos);
        }
    }
    else
    {
        // Get the next line.
        // It is either the description line or "BEGIN TIME STEP".
        readEnsightString(*this, buffer);
    }


    // The buffer string now either contains the description line
    // or "BEGIN TIME STEP"

    if (buffer.starts_with("BEGIN TIME STEP"))
    {
        // Transient single file.
        // File position is now after the "BEGIN TIME STEP" string

        curr_pos = iss.tellg();   // Fallback value

        timeStepFooterBegin_ = getTimeStepFooter(*this, timeStepOffsets_);

        if (timeStepOffsets_.empty())
        {
            // Treat like a single time step
            timeStepOffsets_.resize(1, int64_t(curr_pos));
        }
    }
    else
    {
        // A description line and not "BEGIN TIME STEP"
        // so backtrack to before it was read

        lineNumber(lineNum);        // Restore line number
        iss.seekg(curr_pos);        // Restore file position

        timeStepFooterBegin_ = -1;  // safety
        timeStepOffsets_.clear();   // safety
    }

    DebugInfo<< "Time-steps: " << timeStepOffsets_.size() << endl;

    syncState();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightReadFile::ensightReadFile
(
    const fileName& pathname
)
:
    IFstream(pathname, IOstreamOption::BINARY),  // Start as BINARY
    timeStepFooterBegin_(-1)
{
    init(true);  // detectFormat = true
}


Foam::ensightReadFile::ensightReadFile
(
    const fileName& pathname,
    IOstreamOption::streamFormat fmt
)
:
    IFstream(pathname, fmt),
    timeStepFooterBegin_(-1)
{
    init(false);  // detectFormat = false
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Same as IFstream::readRaw(buf, count)
Foam::Istream& Foam::ensightReadFile::read
(
    char* buf,
    std::streamsize count
)
{
    stdStream().read(buf, count);
    syncState();
    return *this;
}


// TBD
// Foam::Istream& Foam::ensightReadFile::read(word& value)
// {
//     readString(value);
//     string::stripInvalid<word>(value);
//     return *this;
// }


Foam::Istream& Foam::ensightReadFile::read(string& value)
{
    readString(value);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(label& value)
{
    value = getPrimitive<int>(*this);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(float& value)
{
    value = getPrimitive<float>(*this);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(double& value)
{
    value = getPrimitive<float>(*this);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::readKeyword(string& key)
{
    read(key);
    return *this;
}


void Foam::ensightReadFile::readPoints
(
    const label nPoints,
    List<floatVector>& points
)
{
    points.resize_nocopy(nPoints);

    for (auto& p : points)
    {
        read(p.x());
    }
    for (auto& p : points)
    {
        read(p.y());
    }
    for (auto& p : points)
    {
        read(p.z());
    }
}


void Foam::ensightReadFile::readPoints
(
    const label nPoints,
    List<doubleVector>& points
)
{
    points.resize_nocopy(nPoints);

    for (auto& p : points)
    {
        read(p.x());
    }
    for (auto& p : points)
    {
        read(p.y());
    }
    for (auto& p : points)
    {
        read(p.z());
    }
}


bool Foam::ensightReadFile::seekTime(const label timeIndex)
{
    if (timeIndex >= 0 && timeIndex < timeStepOffsets_.size())
    {
        auto& iss = stdStream();

        iss.seekg(timeStepOffsets_[timeIndex]);
        syncState();

        if (debug)
        {
            Info<< "seek time "
                << timeIndex << '/' << nTimes()
                << " offset:" << label(timeStepOffsets_[timeIndex]) << nl;
        }

        return true;
    }

    if (debug)
    {
        Info<< "seek time "
            << timeIndex << '/' << nTimes()
            << " ignored" << nl;
    }

    return false;
}


// ************************************************************************* //
