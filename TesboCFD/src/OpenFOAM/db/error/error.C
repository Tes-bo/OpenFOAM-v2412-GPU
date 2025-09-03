/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
#include "UPstream.H"
#include "StringStream.H"
#include "foamVersion.H"
#include "OSspecific.H"
#include "Enum.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::error::handlerTypes
>
Foam::error::handlerNames
({
    { handlerTypes::DEFAULT, "default" },
    { handlerTypes::IGNORE, "ignore" },
    { handlerTypes::WARN, "warn" },
    { handlerTypes::STRICT, "strict" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::error::master(const label communicator)
{
    // Trap negative value for comm as 'default'. This avoids direct use
    // of UPstream::worldComm which may have not yet been initialised

    return
    (
        UPstream::parRun()
      ? (communicator < 0 ? UPstream::master() : UPstream::master(communicator))
      : true
    );
}


bool Foam::error::warnAboutAge(const int version) noexcept
{
    // No warning for 0 (unversioned) or -ve values (silent versioning)
    return ((version > 0) && (version < foamVersion::api));
}


bool Foam::error::warnAboutAge(const char* what, const int version)
{
    // No warning for 0 (unversioned) or -ve values (silent versioning).
    // Also no warning for (version >= foamVersion::api), which
    // can be used to denote future expiry dates of transition features.

    const bool old = ((version > 0) && (version < foamVersion::api));

    if (old)
    {
        const int months =
        (
            // YYMM -> months
            (12 * (foamVersion::api/100) + (foamVersion::api % 100))
          - (12 * (version/100)  + (version % 100))
        );

        if (version < 1000)
        {
            // For things that predate YYMM versioning (eg, 240 for version 2.4)
            std::cerr
                << "    This " << what << " is very old.\n"
                << std::endl;
        }
        else
        {
            std::cerr
                << "    This " << what << " is " << months << " months old.\n"
                << std::endl;
        }
    }

    return old;
}


bool Foam::error::useAbort()
{
    // FOAM_ABORT env set and contains bool-type value
    return static_cast<bool>(Switch::find(Foam::getEnv("FOAM_ABORT")));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::error::error(const char* title)
:
    std::exception(),
    messageStream(title, messageStream::FATAL),
    functionName_("unknown"),
    sourceFileName_("unknown"),
    sourceFileLineNumber_(0),
    throwing_(false),
    messageStreamPtr_(nullptr)
{}


Foam::error::error(const dictionary& errDict)
:
    std::exception(),
    messageStream(errDict),
    functionName_(errDict.get<string>("functionName")),
    sourceFileName_(errDict.get<string>("sourceFileName")),
    sourceFileLineNumber_(errDict.get<label>("sourceFileLineNumber")),
    throwing_(false),
    messageStreamPtr_(nullptr)
{}


Foam::error::error(const error& err)
:
    std::exception(),
    messageStream(err),
    functionName_(err.functionName_),
    sourceFileName_(err.sourceFileName_),
    sourceFileLineNumber_(err.sourceFileLineNumber_),
    throwing_(err.throwing_),
    messageStreamPtr_(nullptr)
{
    // FIXME: OStringStream copy construct does not adjust tellp and
    // thus the count is wrong! (#3281)
    // // if (err.messageStreamPtr_ && (err.messageStreamPtr_->count() > 0))
    if (err.messageStreamPtr_)
    {
        messageStreamPtr_.reset(new OStringStream(*err.messageStreamPtr_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::error::~error() noexcept
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::error::operator()
(
    const string& functionName
)
{
    functionName_ = functionName;
    sourceFileName_.clear();
    sourceFileLineNumber_ = -1;

    return this->stream();
}


Foam::OSstream& Foam::error::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    functionName_.clear();
    sourceFileName_.clear();

    if (functionName)
    {
        // With nullptr protection
        functionName_.assign(functionName);
    }
    if (sourceFileName)
    {
        // With nullptr protection
        sourceFileName_.assign(sourceFileName);
    }
    sourceFileLineNumber_ = sourceFileLineNumber;

    return this->stream();
}


Foam::OSstream& Foam::error::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


Foam::error::operator Foam::dictionary() const
{
    dictionary errDict;

    string oneLineMessage(message());
    oneLineMessage.replaceAll("\n", " ");

    errDict.add("type", word("Foam::error"));
    errDict.add("message", oneLineMessage);
    errDict.add("function", functionName());
    errDict.add("sourceFile", sourceFileName());
    errDict.add("sourceFileLineNumber", sourceFileLineNumber());

    return errDict;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::error::exiting(const int errNo, const bool isAbort)
{
    if (throwing_)
    {
        if (!isAbort)
        {
            // Make a copy of the error to throw
            error errorException(*this);

            // Reset the message buffer for the next error message
            error::clear();

            throw errorException;
            return;
        }
    }
    else if (JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        JobInfo::shutdown(isAbort || error::useAbort());
    }

    simpleExit(errNo, isAbort);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::error::simpleExit(const int errNo, const bool isAbort)
{
    if (error::useAbort())
    {
        Perr<< nl << *this << nl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        error::printStack(Perr);
        std::abort();
    }
    else if (UPstream::parRun())
    {
        if (isAbort)
        {
            Perr<< nl << *this << nl
                << "\nFOAM parallel run aborting\n" << endl;
            error::printStack(Perr);
            UPstream::abort();
        }
        else
        {
            Perr<< nl << *this << nl
                << "\nFOAM parallel run exiting\n" << endl;
            UPstream::exit(errNo);
        }
    }
    else
    {
        if (isAbort)
        {
            Perr<< nl << *this << nl
                << "\nFOAM aborting\n" << endl;
            error::printStack(Perr);

            #ifdef _WIN32
            std::exit(1);  // Prefer exit() to avoid unnecessary warnings
            #else
            std::abort();
            #endif
        }
        else
        {
            Perr<< nl << *this << nl
                << "\nFOAM exiting\n" << endl;
            std::exit(errNo);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::OSstream& Foam::error::stream()
{
    if (!messageStreamPtr_)
    {
        messageStreamPtr_ = std::make_unique<OStringStream>();
    }
    else if (!messageStreamPtr_->good())
    {
        Perr<< nl
            << "error::stream() : error stream has failed"
            << endl;
        abort();
    }

    return *messageStreamPtr_;
}


Foam::string Foam::error::message() const
{
    if (messageStreamPtr_)
    {
        return messageStreamPtr_->str();
    }

    return string();
}


void Foam::error::clear() const
{
    if (messageStreamPtr_)
    {
        messageStreamPtr_->reset();
    }
}


void Foam::error::exit(const int errNo)
{
    exiting(errNo, false);
}


void Foam::error::abort()
{
    exiting(1, true);
}


void Foam::error::write(Ostream& os, const bool withTitle) const
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

Foam::Ostream& Foam::operator<<(Ostream& os, const error& err)
{
    err.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::error Foam::FatalError("--> FOAM FATAL ERROR: ");


// ************************************************************************* //
