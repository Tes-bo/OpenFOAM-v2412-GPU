/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "codedBase.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "dlLibraryTable.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "Ostream.H"
#include "Time.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedBase, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline void writeEntryIfPresent
(
    Ostream& os,
    const dictionary& dict,
    const word& key
)
{
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);
    if (!eptr)
    {
        // Nothing to do
    }
    else if (eptr->isDict())
    {
        eptr->dict().writeEntry(os);
    }
    else
    {
        const tokenList& toks = eptr->stream();

        if (!toks.empty())  // Could also check that it is a string-type
        {
            os.writeEntry(key, toks[0]);
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::codedBase::writeCodeDict(Ostream& os, const dictionary& dict)
{
    writeEntryIfPresent(os, dict, "codeContext");
    writeEntryIfPresent(os, dict, "codeInclude");
    writeEntryIfPresent(os, dict, "localCode");
    writeEntryIfPresent(os, dict, "code");
    writeEntryIfPresent(os, dict, "codeOptions");
    writeEntryIfPresent(os, dict, "codeLibs");
}


const Foam::dictionary&
Foam::codedBase::codeDict
(
    const objectRegistry& obr,
    const word& dictName
)
{
    auto* dictptr = obr.getObjectPtr<IOdictionary>(dictName);

    if (!dictptr)
    {
        dictptr = new IOdictionary
        (
            IOobject
            (
                dictName,
                obr.time().system(),
                obr,
                IOobject::READ_MODIFIED,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        obr.store(dictptr);
    }

    return *dictptr;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void* Foam::codedBase::loadLibrary
(
    const fileName& libPath,
    const std::string& funcName,
    const dynamicCodeContext& context
) const
{
    // Avoid compilation by loading an existing library

    void* handle = libs().open(libPath, false);

    if (!handle)
    {
        return handle;
    }

    // Verify the loaded version and unload if needed

    // Manual execution of code after loading.
    // This is mandatory for codedBase.

    const bool ok = libs().loadHook(handle, funcName, false);

    if (!ok)
    {
        FatalIOErrorInFunction(context.dict())
            << "Failed symbol lookup " << funcName.c_str() << nl
            << "from " << libPath << nl
            << exit(FatalIOError);

        handle = nullptr;
        if (!libs().close(libPath, false))
        {
            FatalIOErrorInFunction(context.dict())
                << "Failed unloading library " << libPath << nl
                << exit(FatalIOError);
        }
    }

    return handle;
}


void Foam::codedBase::unloadLibrary
(
    const fileName& libPath,
    const std::string& funcName,
    const dynamicCodeContext& context
) const
{
    void* handle = libs().open(libPath, false);

    if (!handle)
    {
        return;
    }

    // Manual execution of code before unloading.
    // This is mandatory for codedBase.

    const bool ok = libs().unloadHook(handle, funcName, false);

    if (!ok)
    {
        IOWarningInFunction(context.dict())
            << "Failed looking up symbol " << funcName << nl
            << "from " << libPath << nl;
    }

    if (!libs().close(libPath, false))
    {
        FatalIOErrorInFunction(context.dict())
            << "Failed unloading library " << libPath << nl
            << exit(FatalIOError);
    }
}


void Foam::codedBase::createLibrary
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // No library?
    // The overall master should compile it (so only it needs a compiler) and
    // - see if it appears at the other processors
    // - or copy it across if it does not

    // Indicates NFS filesystem
    const bool isNFS = (IOobject::fileModificationSkew > 0);

    if
    (
        (UPstream::master() || !isNFS)
     && !dynCode.upToDate(context)
    )
    {
        // Filter with this context
        dynCode.reset(context);

        this->prepare(dynCode, context);

        if (!dynCode.copyOrCreateFiles(true))
        {
            FatalIOErrorInFunction(context.dict())
                << "Failed writing files for" << nl
                << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }

        if (!dynCode.wmakeLibso())
        {
            FatalIOErrorInFunction(context.dict())
                << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }

    if (isNFS)
    {
        // Wait for compile to finish before attempting filesystem polling
        UPstream::barrier(UPstream::worldComm);
    }

    const fileName libPath = dynCode.libPath();

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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::codedBase::setCodeContext(const dictionary& dict)
{
    context_.setCodeContext(dict);
}


void Foam::codedBase::append(const std::string& str)
{
    context_.append(str);
}


void Foam::codedBase::updateLibrary
(
    const word& name,
    const dynamicCodeContext& context
) const
{
    dynamicCode::checkSecurity
    (
        "codedBase::updateLibrary()",
        context.dict()
    );

    // codeName: name + _<sha1>
    // codeDir : name
    dynamicCode dynCode
    (
        name + context.sha1().str(true),
        name
    );

    const fileName libPath = dynCode.libPath();


    // The correct library was already loaded => we are done
    if (libs().findLibrary(libPath))
    {
        return;
    }

    DetailInfo
        << "Using dynamicCode for " << this->description().c_str()
        << " at line " << context.dict().startLineNumber()
        << " in " << context.dict().name() << endl;


    // Remove instantiation of fvPatchField provided by library
    this->clearRedirect();

    // May need to unload old library
    unloadLibrary(oldLibPath_, dlLibraryTable::basename(oldLibPath_), context);

    // Try loading an existing library (avoid compilation when possible)
    void* lib = loadLibrary(libPath, dynCode.codeName(), context);

    if (returnReduceOr(lib == nullptr))
    {
        if (lib)
        {
            // Ensure consistency
            unloadLibrary(libPath, dlLibraryTable::basename(libPath), context);
        }

        createLibrary(dynCode, context);

        loadLibrary(libPath, dynCode.codeName(), context);
    }

    // Retain for future reference
    oldLibPath_ = libPath;
}


void Foam::codedBase::updateLibrary
(
    const word& name,
    const dictionary& dict
) const
{
    updateLibrary(name, dynamicCodeContext(dict));
}


void Foam::codedBase::updateLibrary(const word& name) const
{
    if (context_.good())
    {
        updateLibrary(name, context_);
    }
    else
    {
        updateLibrary(name, dynamicCodeContext(this->codeDict()));
    }
}


// ************************************************************************* //
