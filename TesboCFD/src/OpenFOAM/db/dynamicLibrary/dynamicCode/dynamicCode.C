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

#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "dlLibraryTable.H"
#include "argList.H"
#include "stringOps.H"
#include "Fstream.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "etcFiles.H"
#include "dictionary.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::dynamicCode::allowSystemOperations
(
    Foam::debug::infoSwitch("allowSystemOperations", 0)
);


const Foam::word Foam::dynamicCode::codeTemplateEnvName
    = "FOAM_CODE_TEMPLATES";

const Foam::fileName Foam::dynamicCode::codeTemplateDirName
    = "codeTemplates/dynamicCode";

const char* const Foam::dynamicCode::targetLibDir
    = "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib";

const char* const Foam::dynamicCode::topDirName
    = "dynamicCode";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::dynamicCode::checkSecurity
(
    const char* title,
    const dictionary& dict
)
{
    if (isAdministrator())
    {
        FatalIOErrorInFunction(dict)
            << "This code should not be executed by someone"
            << " with administrator rights for security reasons." << nl
            << "It generates a shared library which is loaded using dlopen"
            << nl << endl
            << exit(FatalIOError);
    }

    if (!allowSystemOperations)
    {
        FatalIOErrorInFunction(dict)
            << "Loading shared libraries using case-supplied code may have"
            << " been disabled" << nl
            << "by default for security reasons." << nl
            << "If you trust the code, you may enable this by adding"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "to the InfoSwitches setting in the system controlDict." << nl
            << "The system controlDict is any of" << nl << nl
            << "    ~/.OpenFOAM/" << foamVersion::api << "/controlDict" << nl
            << "    ~/.OpenFOAM/controlDict" << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl << endl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicCode::copyAndFilter
(
    ISstream& is,
    OSstream& os,
    const HashTable<string>& mapping
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "Failed opening for reading " << is.name()
            << exit(FatalError);
    }

    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed writing " << os.name()
            << exit(FatalError);
    }

    // Copy file while rewriting $VARS and ${VARS}
    string line;
    do
    {
        is.getLine(line);

        // Expand according to HashTable mapping, not the environment.
        // Expanding according to env variables might cause too many
        // surprises
        stringOps::inplaceExpand(line, mapping);
        os.writeQuoted(line, false) << nl;
    }
    while (is.good());
}


bool Foam::dynamicCode::resolveTemplates
(
    const UList<fileName>& templateNames,
    DynamicList<fileName>& resolvedFiles,
    DynamicList<fileName>& badFiles
)
{
    // Try to get template from FOAM_CODE_TEMPLATES
    const fileName templateDir(Foam::getEnv(codeTemplateEnvName));

    bool allOkay = true;
    for (const fileName& templateName : templateNames)
    {
        fileName file;
        if (!templateDir.empty() && isDir(templateDir))
        {
            file = templateDir/templateName;
            if (!isFile(file, false))
            {
                file.clear();
            }
        }

        // Not found - fallback to <etc> expansion
        if (file.empty())
        {
            file = findEtcFile(codeTemplateDirName/templateName);
        }

        if (file.empty())
        {
            badFiles.push_back(templateName);
            allOkay = false;
        }
        else
        {
            resolvedFiles.push_back(file);
        }
    }

    return allOkay;
}


bool Foam::dynamicCode::writeCommentSHA1(Ostream& os) const
{
    const auto iter = filterVars_.cfind("SHA1sum");

    if (iter.good())
    {
        os  << "/* dynamicCode:\n * SHA1 = ";
        os.writeQuoted(iter.val(), false) << "\n */\n";
    }

    return iter.good();
}


bool Foam::dynamicCode::createMakeFiles() const
{
    // Create Make/files
    if (compileFiles_.empty())
    {
        return false;
    }

    const fileName dstFile(this->codePath()/"Make/files");

    // Create dir
    mkDir(dstFile.path());

    OFstream os(dstFile);
    //Debug: Info<< "Writing to " << dstFile << endl;
    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed writing " << dstFile
            << exit(FatalError);
    }

    writeCommentSHA1(os);

    // Write compile files
    for (const fileName& file : compileFiles_)
    {
        os.writeQuoted(file, false) << nl;
    }

    os  << nl
        << targetLibDir
        << "/lib" << codeName_.c_str() << nl;

    return true;
}


bool Foam::dynamicCode::createMakeOptions() const
{
    // Create Make/options
    if (compileFiles_.empty() || makeOptions_.empty())
    {
        return false;
    }

    const fileName dstFile(this->codePath()/"Make/options");

    // Create dir
    mkDir(dstFile.path());

    OFstream os(dstFile);
    //Debug: Info<< "Writing to " << dstFile << endl;
    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed writing " << dstFile
            << exit(FatalError);
    }

    writeCommentSHA1(os);
    os.writeQuoted(makeOptions_, false) << nl;

    return true;
}


bool Foam::dynamicCode::writeDigest(const SHA1Digest& sha1) const
{
    const fileName file = digestFile();
    mkDir(file.path());

    OFstream os(file);
    sha1.write(os, true) << nl;

    return os.good();
}


bool Foam::dynamicCode::writeDigest(const std::string& sha1) const
{
    const fileName file = digestFile();
    mkDir(file.path());

    OFstream os(file);
    os  << '_';
    os.writeQuoted(sha1, false) << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCode::dynamicCode(const word& codeName, const word& codeDirName)
:
    codeRoot_(argList::envGlobalPath()/topDirName),
    libSubDir_(stringOps::expand("platforms/${WM_OPTIONS}/lib")),
    codeName_(codeName),
    codeDirName_(codeDirName)
{
    if (codeDirName_.empty())
    {
        codeDirName_ = codeName_;
    }

    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::dynamicCode::codeRelPath() const
{
    return topDirName/codeDirName_;
}


Foam::fileName Foam::dynamicCode::libPath() const
{
    return codeRoot_/libSubDir_/dlLibraryTable::fullname(codeName_);
}


Foam::fileName Foam::dynamicCode::libRelPath() const
{
    return codeRelPath()/libSubDir_/dlLibraryTable::fullname(codeName_);
}


void Foam::dynamicCode::clear()
{
    compileFiles_.clear();
    copyFiles_.clear();
    createFiles_.clear();
    filterVars_.clear();
    filterVars_.set("typeName", codeName_);
    filterVars_.set("SHA1sum", SHA1Digest().str());

    // Default Make/options
    makeOptions_ =
        "EXE_INC = -g\n"
        "\n\nLIB_LIBS = ";
}


void Foam::dynamicCode::reset
(
    const dynamicCodeContext& context
)
{
    clear();
    setFilterContext(context);
}


void Foam::dynamicCode::addCompileFile(const fileName& name)
{
    compileFiles_.push_back(name);
}


void Foam::dynamicCode::addCopyFile(const fileName& name)
{
    copyFiles_.push_back(name);
}


void Foam::dynamicCode::addCreateFile
(
    const fileName& name,
    const std::string& fileContents
)
{
    createFiles_.emplace_back(name, fileContents);
}


void Foam::dynamicCode::setFilterContext
(
    const dynamicCodeContext& context
)
{
    filterVars_.set("localCode", context.localCode());
    filterVars_.set("code", context.code());
    filterVars_.set("codeInclude", context.include());
    filterVars_.set("SHA1sum", context.sha1().str());
}


void Foam::dynamicCode::setFilterVariable
(
    const word& key,
    const std::string& value
)
{
    filterVars_.set(key, value);
}


void Foam::dynamicCode::setMakeOptions(const std::string& content)
{
    makeOptions_ = content;
}


bool Foam::dynamicCode::copyOrCreateFiles(const bool verbose) const
{
    if (verbose)
    {
        DetailInfo
            << "Creating new library in " << this->libRelPath() << endl;
    }

    const label nFiles = compileFiles_.size() + copyFiles_.size();

    DynamicList<fileName> resolvedFiles(nFiles);
    DynamicList<fileName> badFiles(nFiles);

    // Resolve template, or add to bad-files
    resolveTemplates(compileFiles_, resolvedFiles, badFiles);
    resolveTemplates(copyFiles_, resolvedFiles, badFiles);

    if (!badFiles.empty())
    {
        FatalErrorInFunction
            << "Could not find code template(s): "
            << badFiles << nl
            << "Under the $" << codeTemplateEnvName
            << " directory or via the <etc>/"
            << codeTemplateDirName << " expansion"
            << exit(FatalError);
    }



    // Create dir
    const fileName outputDir = this->codePath();

    // Create dir
    mkDir(outputDir);

    // Copy/filter files
    for (const fileName& srcFile : resolvedFiles)
    {
        const fileName dstFile(outputDir/srcFile.name());

        IFstream is(srcFile);
        //Debug: Info<< "Reading from " << is.name() << endl;
        if (!is.good())
        {
            FatalErrorInFunction
                << "Failed opening " << srcFile
                << exit(FatalError);
        }

        OFstream os(dstFile);
        //Debug: Info<< "Writing to " << dstFile.name() << endl;
        if (!os.good())
        {
            FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
        }

        // Copy lines while expanding variables
        copyAndFilter(is, os, filterVars_);
    }


    // Create files:
    for (const auto& content : createFiles_)
    {
        const fileName dstFile(outputDir/stringOps::expand(content.first()));

        mkDir(dstFile.path());
        OFstream os(dstFile);
        //Debug: Info<< "Writing to " << content.first() << endl;
        if (!os.good())
        {
            FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
        }
        os.writeQuoted(content.second(), false) << nl;
    }


    // Create Make/files + Make/options
    createMakeFiles();
    createMakeOptions();

    writeDigest(filterVars_["SHA1sum"]);

    return true;
}


bool Foam::dynamicCode::wmakeLibso() const
{
    stringList cmd({"wmake", "-s", "libso", this->codePath()});

    // NOTE: could also resolve wmake command explicitly
    //   cmd[0] = stringOps::expand("$WM_PROJECT_DIR/wmake/wmake");

    // This can take a bit longer, so report that we are starting wmake
    // Even with details turned off, we want some feedback

    OSstream& os = (Foam::infoDetailLevel > 0 ? Info : InfoErr);
    os  << "Invoking wmake libso " << this->codePath().c_str() << endl;

    if (Foam::system(cmd) == 0)
    {
        return true;
    }

    return false;
}


bool Foam::dynamicCode::upToDate(const SHA1Digest& sha1) const
{
    const fileName file = digestFile();

    if (!exists(file, false) || SHA1Digest(IFstream(file)()) != sha1)
    {
        return false;
    }

    return true;
}


bool Foam::dynamicCode::upToDate(const dynamicCodeContext& context) const
{
    return upToDate(context.sha1());
}


// * * * * * * * * * * * * * * * Synchronisation * * * * * * * * * * * * * * //

void Foam::dynamicCode::waitForFile
(
    const fileName& file,
    const dictionary& contextDict
)
{
    const int debug = 0;

    if (!UPstream::parRun())
    {
        return;
    }

    // If library has just been compiled on the master, the other nodes
    // need to pick this library up through NFS, which likely has delays
    // in it.

    // We do this by just polling a few times using the
    // fileModificationSkew.

    off_t localSize = Foam::fileSize(file);
    off_t masterSize = localSize;
    Pstream::broadcast(masterSize);

    for
    (
        label iter = 0;
        (
            iter < IOobject::maxFileModificationPolls
         && IOobject::fileModificationSkew > 0
        );
        ++iter
    )
    {
        DebugPout
            << "Processor " << UPstream::myProcNo()
            << " masterSize:" << masterSize
            << " localSize:" << localSize << endl;

        if (localSize == masterSize)
        {
            return;
        }
        if (localSize > masterSize)
        {
            FatalIOErrorInFunction(contextDict)
                << "Excessive size when reading (NFS mounted) library "
                << nl << file << nl
                << "on processor " << UPstream::myProcNo()
                << " detected size " << localSize
                << " whereas master size is " << masterSize
                << " bytes." << nl
                << "If your case is NFS mounted increase"
                << " fileModificationSkew or maxFileModificationPolls;"
                << nl << "If your case is not NFS mounted"
                << " (so distributed) set fileModificationSkew"
                << " to 0"
                << exit(FatalIOError);
        }
        else
        {
            DebugPout
                << "Local file " << file
                << " not of same size (" << localSize
                << ") as master ("
                << masterSize << "). Waiting for "
                << IOobject::fileModificationSkew
                << " seconds." << endl;

            Foam::sleep(IOobject::fileModificationSkew);

            // Recheck local size
            localSize = Foam::fileSize(file);
        }
    }

    // Finished doing iterations. Do final check
    if (localSize != masterSize)
    {
        FatalIOErrorInFunction(contextDict)
            << "Cannot read (NFS mounted) library:" << nl
            << file << nl
            << "on processor " << UPstream::myProcNo()
            << " detected size " << localSize
            << " whereas master size is " << masterSize
            << " bytes." << nl
            << "If your case is NFS mounted increase"
            << " fileModificationSkew or maxFileModificationPolls;" << nl
            << "If your case is not NFS mounted"
            << " (so distributed) set fileModificationSkew"
            << " to 0"
            << exit(FatalIOError);
    }

    DebugPout
        << "Processor " << UPstream::myProcNo()
        << " masterSize:" << masterSize
        << " localSize:" << localSize
        << " ... after waiting" << endl;
}


// ************************************************************************* //
