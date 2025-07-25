/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::messageStream

Description
    Handle output messages in a simple, consistent stream-based manner.

    The messageStream class is globally instantiated with a title
    string and a severity (which controls the program termination),
    optionally with a max number of errors before termination.

    Errors, messages and other data are sent to the messageStream class in
    the standard manner.

    For parallel applications, the output for 'standard' messages
    (Info, Warnings) is effectively suppressed on all sub-processes,
    which results in a single output message instead of a flood of output
    messages from each process.  The error type of messages do, however,
    retain output on all processes, which ensures that parallel termination
    occurs correctly and the source of the problem is properly traceable to
    the originating processor.

SourceFiles
    messageStream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_messageStream_H
#define Foam_messageStream_H

#include "label.H"
#include "word.H"
#include <string>
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOstream;
class OSstream;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class messageStream Declaration
\*---------------------------------------------------------------------------*/

class messageStream
{
public:

    //- Message type, error severity flags
    enum errorSeverity : int
    {
        // Serial-only output:
        INFO = 1,       //!< General information output (stdout)
        WARNING,        //!< Warning of possible problem.

        // Parallel-aware output:
        SERIOUS,        //!< A serious problem - eg, data corruption.
        FATAL,          //!< A fatal error.

        //! Bitmask for stderr output (for the above enums)
        USE_STDERR = 0x80
    };


protected:

    // Protected Data

        //- The title of this error type
        string title_;

        //- The message type / error severity, possibly with USE_STDERR mask
        int severity_;

        //- The maximum number of errors before program termination
        int maxErrors_;

        //- The current number of errors counted
        int errorCount_;


public:

    // Static Data

        //- The output level (verbosity) of messages
        //
        //  - level == 0 : suppress all output
        //  - level == 1 : normal output
        //  - level >= 2 : report source file name and line number if available
        //
        // \note The default level is normally 2.
        static int level;

        //- The output redirection of messages
        //
        //  - redirect == 2 : use stderr instead of stdout
        static int redirect;


    // Constructors

        //- Construct from components
        messageStream
        (
            const char* title,
            errorSeverity severity,
            int maxErrors = 0,
            bool use_stderr = false
        );

        //- Construct from components
        messageStream
        (
            const std::string& title,
            errorSeverity severity,
            int maxErrors = 0,
            bool use_stderr = false
        )
        :
            messageStream(title.c_str(), severity, maxErrors, use_stderr)
        {}


        //- Construct from dictionary as Fatal, extracting 'title'.
        explicit messageStream(const dictionary& dict);


    // Member Functions

        //- The title of this error type
        const string& title() const noexcept
        {
            return title_;
        }

        //- The maximum number of errors before program termination
        int maxErrors() const noexcept
        {
            return maxErrors_;
        }

        //- Set the maximum number of errors before program termination
        //  \return the previous value for maxErrors
        int maxErrors(int nErrors) noexcept
        {
            int old = maxErrors_;
            maxErrors_ = nErrors;
            return old;
        }


    // Output

        //- Return OSstream for output operations.
        OSstream& stream
        (
            //! An alternative output stream (serial-only)
            OSstream* alternative = nullptr
        );

        //- Return OSstream for output operations on the master process only,
        //- Snull on other processes.
        OSstream& masterStream(const label communicator);

        //- Return std::ostream for output operations.
        std::ostream& stdStream();

        //- Report deprecation (after specified API version) with
        //- 'From function-name, source file, line number'.
        //  \return OSstream for further operations
        OSstream& deprecated
        (
            const int afterVersion,
            const char* functionName = nullptr,
            const char* sourceFileName = nullptr,
            const int sourceFileLineNumber = 0
        );


        //- Implicit cast to OSstream for << operations
        operator OSstream&()
        {
            return this->stream();
        }

        //- Explicitly convert to OSstream for << operations
        OSstream& operator()()
        {
            return this->stream();
        }

        //- Report 'From function-name'
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const std::string& functionName
        );

        //- Report 'From function-name, source file, line number'
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Report 'From function-name, source file, line number'
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const std::string& functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Report 'From function-name, source file, line number'
        //- as well as io-file name and location
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const std::string& ioFileName,
            const label ioStartLineNumber = -1,
            const label ioEndLineNumber = -1
        );

        //- Report 'From function-name, source file, line number'
        //- as well as io-file name and location
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream&  //!< Provides name and lineNumber
        );

        //- Report 'From function-name, source file, line number'
        //- as well as io-file name and location
        //  \return OSstream for further operations
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const dictionary&  //!< Provides name, startLine, endLine
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error declarations: defined in messageStream.C

//- Global for selective suppression of Info output.
//  This is normally accessed implicitly via the DetailInfo macro and is often
//  associated with applications with suppressed banners. For example,
//
//  \code
//      DetailInfo << "Hello, I'm running from program xyz" << nl;
//      Info<< "Found ... invalid items" << nl;
//  \endcode
//
//  The values are normally 0 or a positive value.
//  \note This flag is initialized to 1 by default.
extern int infoDetailLevel;

//- Information stream (stdout output on master, null elsewhere)
extern messageStream Info;

//- Information stream (stderr output on master, null elsewhere)
extern messageStream InfoErr;

//- Warning stream (stdout output on master, null elsewhere),
//- with additional 'FOAM Warning' header text.
extern messageStream Warning;

//- Error stream (stdout output on all processes),
//- with additional 'FOAM Serious Error' header text.
extern messageStream SeriousError;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Convenience macros to add the file name and line number to the function name

// Compiler provided function name string:
//     for gcc-compatible compilers use __PRETTY_FUNCTION__
//     otherwise use the standard __func__
#ifdef __GNUC__
    #define FUNCTION_NAME __PRETTY_FUNCTION__
#else
    #define FUNCTION_NAME __func__
#endif


//- Report an error message using Foam::SeriousError
//  for functionName in file __FILE__ at line __LINE__
#define SeriousErrorIn(functionName)                                           \
    ::Foam::SeriousError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam::SeriousError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define SeriousErrorInFunction SeriousErrorIn(FUNCTION_NAME)


//- Report an IO error message using Foam::SeriousError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SeriousIOErrorIn(functionName, ios)                                    \
    ::Foam::SeriousError((functionName), __FILE__, __LINE__, ios)

//- Report an IO error message using Foam::SeriousError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SeriousIOErrorInFunction(ios) SeriousIOErrorIn(FUNCTION_NAME, ios)


//- Report a warning using Foam::Warning
//  for functionName in file __FILE__ at line __LINE__
#define WarningIn(functionName)                                                \
    ::Foam::Warning((functionName), __FILE__, __LINE__)

//- Report a warning using Foam::Warning
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define WarningInFunction WarningIn(FUNCTION_NAME)

//- Report a warning using Foam::Warning
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define DeprecatedInFunction(afterVersion)                                     \
    ::Foam::Warning.deprecated(afterVersion, FUNCTION_NAME, __FILE__, __LINE__)


//- Report an IO warning using Foam::Warning
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOWarningIn(functionName, ios)                                         \
    ::Foam::Warning((functionName), __FILE__, __LINE__, (ios))

//- Report an IO warning using Foam::Warning
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOWarningInFunction(ios) IOWarningIn(FUNCTION_NAME, ios)


//- Report an information message using Foam::Info
//  for functionName in file __FILE__ at line __LINE__
#define InfoIn(functionName)                                                   \
    ::Foam::Info((functionName), __FILE__, __LINE__)

//- Report an information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define InfoInFunction InfoIn(FUNCTION_NAME)

//- Report using Foam::Pout with functionName: prefix
#define PoutIn(functionName)                                                   \
    ::Foam::Pout << (functionName) << ':'

//- Report using Foam::Pout with FUNCTION_NAME prefix
#define PoutInFunction PoutIn(FUNCTION_NAME)

//- Write to Foam::Info if the Foam::infoDetailLevel is +ve non-zero (default)
#define DetailInfo                                                             \
    if (::Foam::infoDetailLevel > 0) Info

//- Report write to Foam::Info if the local log switch is true
#define Log                                                                    \
    if (log) Info

//- Report write to Foam::Info if the class log switch is true
#define Log_                                                                   \
    if (this->log) Info


//- Report an IO information message using Foam::Info
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOInfoIn(functionName, ios)                                            \
    ::Foam::Info((functionName), __FILE__, __LINE__, (ios))

//- Report an IO information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOInfoInFunction(ios) IOInfoIn(FUNCTION_NAME, ios)


//- Report an information message using Foam::Info
//  if the local debug switch is true
#define DebugInfo                                                              \
    if (debug) Info

//- Report an information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  if the local debug switch is true
#define DebugInFunction                                                        \
    if (debug) InfoInFunction

//- Report an information message using Foam::Pout
//  if the local debug switch is true
#define DebugPout                                                              \
    if (debug) Pout

//- Report an information message using Foam::Pout
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  if the local debug switch is true
#define DebugPoutInFunction                                                    \
    if (debug) PoutInFunction

//- Report a variable name and value
//  using Foam::Pout in file __FILE__ at line __LINE__
#define DebugVar(var)                                                          \
{                                                                              \
    ::Foam::string oldPrefix(::Foam::Pout.prefix());                           \
    ::Foam::Pout<< "["<< __FILE__ << ":" << __LINE__ << "] ";                  \
    ::Foam::Pout.prefix() = oldPrefix + #var " ";                              \
    ::Foam::Pout<< var << ::Foam::endl;                                        \
    ::Foam::Pout.prefix() = oldPrefix;                                         \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
