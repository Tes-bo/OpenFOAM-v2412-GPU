/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "codeStream.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "StringStream.H"
#include "Time.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        dictionaryIstream,
        codeStream
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream,
        codeStream
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable& Foam::functionEntries::codeStream::libs
(
    const dictionary& dict
)
{
    return static_cast<const baseIOdictionary&>(dict.topDict()).time().libs();
}


bool Foam::functionEntries::codeStream::doingMasterOnlyReading
(
    const dictionary& dict
)
{
    // Fallback value
    bool masterOnly = regIOobject::masterOnlyReading;

    const auto* rioPtr = isA<regIOobject>(dict.topDict());

    if (rioPtr)
    {
        masterOnly = rioPtr->global();
    }

    DebugPout
        << "codeStream : " << (rioPtr ? "IO" : "plain")
        << " dictionary:" << dict.name()
        << " master-only-reading:" << masterOnly << endl;

    return masterOnly;
}


Foam::functionEntries::codeStream::streamingFunctionType
Foam::functionEntries::codeStream::getFunction
(
    const dictionary& parentDict,
    const dictionary& codeDict
)
{
    // get code, codeInclude, codeOptions
    dynamicCodeContext context(codeDict);

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    std::string sha1Str(context.sha1().str(true));
    dynamicCode dynCode("codeStream" + sha1Str, sha1Str);


    const dictionary& topDict = parentDict.topDict();
    const bool masterOnly = doingMasterOnlyReading(topDict);

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    const fileName libPath = dynCode.libPath();

    // See if library is loaded
    void* lib = nullptr;

    if (isA<baseIOdictionary>(topDict))
    {
        lib = libs(parentDict).findLibrary(libPath);
    }

    // nothing loaded
    // avoid compilation if possible by loading an existing library
    if (!lib)
    {
        DetailInfo
            << "Using #codeStream with " << libPath << endl;

        if (isA<baseIOdictionary>(topDict))
        {
            // Cached access to libs, with cleanup upon termination
            lib = libs(parentDict).open(libPath, false);
        }
        else
        {
            // Uncached opening of libPath. Do not complain if cannot be loaded
            lib = Foam::dlOpen(libPath, false);
        }
    }


    // Indicates NFS filesystem
    const bool isNFS = (IOobject::fileModificationSkew > 0);

    // Create library if required
    if
    (
        lib == nullptr
     && (UPstream::master() || !isNFS)
     && !dynCode.upToDate(context)
    )
    {
        // Filter with this context
        dynCode.reset(context);

        // Compile filtered C template
        dynCode.addCompileFile(codeTemplateC);

        // define Make/options
        dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
          + context.options()
          + "\n\nLIB_LIBS = \\\n"
            "    -lOpenFOAM \\\n"
          + context.libs()
        );

        if (!dynCode.copyOrCreateFiles(true))
        {
            FatalIOErrorInFunction(parentDict)
                << "Failed writing files for" << nl
                << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }

        if (!dynCode.wmakeLibso())
        {
            FatalIOErrorInFunction(parentDict)
                << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    //- Only block if we're not doing master-only reading.
    //  (flag set by regIOobject::read, baseIOdictionary constructor)
    if (!masterOnly && returnReduceOr(lib == nullptr))
    {
        // Broadcast to distributed masters
        if (fileHandler().distributed())
        {
            fileHandler().broadcastCopy
            (
                UPstream::worldComm,
                UPstream::master(fileHandler().comm()),
                libPath,
                libPath
            );
        }

        dynamicCode::waitForFile(libPath, context.dict());
    }

    if (!lib)
    {
        if (isA<baseIOdictionary>(topDict))
        {
            lib = libs(parentDict).open(libPath, false);
        }
        else
        {
            lib = Foam::dlOpen(libPath, false);
        }
    }


    if (masterOnly ? !lib : returnReduceOr(!lib))
    {
        FatalIOErrorInFunction(parentDict)
            << "Failed loading library " << dynCode.libRelPath()
            << " on some processors."
            << "Did you add all libraries to the 'libs' entry"
            << " in system/controlDict?"
            << exit(FatalIOError);
    }


    // Find the function handle in the library
    streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            Foam::dlSym(lib, dynCode.codeName())
        );


    if (!function)
    {
        FatalIOErrorInFunction(parentDict)
            << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << dynCode.libRelPath()
            << exit(FatalIOError);
    }

    return function;
}


Foam::string Foam::functionEntries::codeStream::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    DetailInfo
        << "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.relativeName() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::evaluate(..)",
        parentDict
    );

    // Get code dictionary
    dictionary codeDict("#codeStream", parentDict, is);

    // Use function to write stream
    OStringStream os(is.format());

    streamingFunctionType function = getFunction(parentDict, codeDict);
    (*function)(os, parentDict);

    // Return evaluated content as string
    return os.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    entry.read(parentDict, result);

    return true;
}


bool Foam::functionEntries::codeStream::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    parentDict.read(result);

    return true;
}


// ************************************************************************* //
