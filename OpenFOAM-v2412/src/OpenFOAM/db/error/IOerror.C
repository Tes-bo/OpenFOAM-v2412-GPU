/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "fileName.H"
#include "dictionary.H"
#include "JobInfo.H"
#include "Pstream.H"
#include "StringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOerror::IOerror(const char* title)
:
    error(title),
    ioFileName_("unknown"),
    ioStartLineNumber_(-1),
    ioEndLineNumber_(-1)
{}


Foam::IOerror::IOerror(const dictionary& errDict)
:
    error(errDict),
    ioFileName_(errDict.get<string>("ioFileName")),
    ioStartLineNumber_(errDict.get<label>("ioStartLineNumber")),
    ioEndLineNumber_(errDict.get<label>("ioEndLineNumber"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOerror::~IOerror() noexcept
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const string& ioFileName,
    const label ioStartLineNumber,
    const label ioEndLineNumber
)
{
    OSstream& os = error::operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber
    );

    ioFileName_ = ioFileName;
    ioStartLineNumber_ = ioStartLineNumber;
    ioEndLineNumber_ = ioEndLineNumber;

    return os;
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        ioStream.relativeName(),
        ioStream.lineNumber(),
        -1  // No known endLineNumber
    );
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const dictionary& dict
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        dict.relativeName(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


Foam::OSstream& Foam::IOerror::operator()
(
    const std::string& where,
    const IOstream& ioStream
)
{
    return operator()
    (
        where.c_str(),
        "",     // No source file
        -1,     // Non-zero to ensure 'where' is reported
        ioStream.relativeName(),
        ioStream.lineNumber(),
        -1      // No known endLineNumber
    );
}


Foam::OSstream& Foam::IOerror::operator()
(
    const std::string& where,
    const dictionary& dict
)
{
    return operator()
    (
        where.c_str(),
        "",     // No source file
        -1,     // Non-zero to ensure 'where' is reported
        dict.relativeName(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


void Foam::IOerror::SafeFatalIOError
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream,
    const std::string& msg
)
{
    if (JobInfo::constructed)
    {
        FatalIOError
        (
            functionName,
            sourceFileName,
            sourceFileLineNumber,
            ioStream
        )   << msg.c_str() << Foam::exit(FatalIOError);
    }
    else
    {
        // Without (openfoam=API patch=NN) since it is rarely used
        std::cerr
            << nl
            << "--> FOAM FATAL IO ERROR:" << nl
            << msg << nl
            << "file: " << ioStream.relativeName()
            << " at line " << ioStream.lineNumber() << '.' << nl << nl
            << "    From " << functionName << nl
            << "    in file " << sourceFileName
            << " at line " << sourceFileLineNumber << '.' << std::endl;
        std::exit(1);
    }
}


Foam::IOerror::operator Foam::dictionary() const
{
    dictionary errDict(error::operator dictionary());

    errDict.add("type", word("Foam::IOerror"), true);  // overwrite
    errDict.add("ioFileName", ioFileName());
    errDict.add("ioStartLineNumber", ioStartLineNumber());
    errDict.add("ioEndLineNumber", ioEndLineNumber());

    return errDict;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IOerror::exiting(const int errNo, const bool isAbort)
{
    if (throwing_)
    {
        if (!isAbort)
        {
            // Make a copy of the error to throw
            IOerror errorException(*this);

            // Reset the message buffer for the next error message
            error::clear();

            throw errorException;
            return;
        }
    }
    else if (JobInfo::constructed)
    {
        jobInfo.add("FatalIOError", operator dictionary());
        JobInfo::shutdown(isAbort || error::useAbort());
    }

    simpleExit(errNo, isAbort);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IOerror::exit(const int)
{
    exiting(1, false);
}


void Foam::IOerror::abort()
{
    exiting(1, true);
}


void Foam::IOerror::write(Ostream& os, const bool withTitle) const
{
    if (os.bad())
    {
        return;
    }

    os  << nl;
    if (withTitle && !title().empty())
    {
        os  << title().c_str()
            << "(openfoam-" << foamVersion::api;

        if (foamVersion::patched())
        {
            // Patch-level, when defined
            os  << " patch=" << foamVersion::patch.c_str();
        }
        os  << ')' << nl;
    }
    os  << message().c_str();


    if (!ioFileName().empty())
    {
        os  << nl << nl
            << "file: " << ioFileName().c_str();

        if (ioStartLineNumber() >= 0)
        {
            os  << " at line " << ioStartLineNumber();
            if (ioStartLineNumber() < ioEndLineNumber())
            {
                os  << " to " << ioEndLineNumber();
            }
            os  << '.';
        }
    }


    const label lineNo = sourceFileLineNumber();

    if (messageStream::level >= 2 && lineNo && !functionName().empty())
    {
        os  << nl << nl
            << "    From " << functionName().c_str() << nl;

        if (!sourceFileName().empty())
        {
            os << "    in file " << sourceFileName().c_str();

            if (lineNo > 0)
            {
                os  << " at line " << lineNo << '.';
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOerror& err)
{
    err.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::IOerror Foam::FatalIOError("--> FOAM FATAL IO ERROR: ");


// ************************************************************************* //
