/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "functionObjectList.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "profiling.H"
#include "argList.H"
#include "timeControlFunctionObject.H"
#include "dictionaryEntry.H"
#include "stringOps.H"
#include "Switch.H"
#include "Tuple2.H"
#include "etcFiles.H"
#include "IOdictionary.H"
#include "Pstream.H"
#include "OSspecific.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

//- Max number of warnings (per functionObject)
static constexpr const unsigned maxWarnings = 10u;

Foam::fileName Foam::functionObjectList::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Mimic exit handling of the error class
    static void exitNow(const error& err)
    {
        if (error::useAbort())
        {
            Perr<< nl << err << nl
                << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
            error::printStack(Perr);
            std::abort();
        }
        else if (UPstream::parRun())
        {
            Perr<< nl << err << nl
                << "\nFOAM parallel run exiting\n" << endl;
            UPstream::exit(1);
        }
        else
        {
            Perr<< nl << err << nl
                << "\nFOAM exiting\n" << endl;
            std::exit(1);
        }
    }

} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::functionObjectList::createPropertiesDict() const
{
    // Cannot set the properties dictionary on construction since Time has not
    // been fully initialised
    propsDictPtr_.reset
    (
        new functionObjects::properties
        (
            IOobject
            (
                "functionObjectProperties",
                time_.timeName(),
                "uniform"/word("functionObjects"),
                time_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        )
    );
}


void Foam::functionObjectList::createOutputRegistry() const
{
    objectsRegistryPtr_.reset
    (
        new objectRegistry
        (
            IOobject
            (
                "functionObjectObjects",
                time_.timeName(),
                time_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        )
    );
}


Foam::autoPtr<Foam::functionObject> Foam::functionObjectList::remove
(
    const word& key,
    label& oldIndex
)
{
    autoPtr<functionObject> oldptr;

    auto iter = indices_.find(key);  // Index of existing functionObject

    if (iter.good())
    {
        oldIndex = *iter;

        // Remove pointer from the old list
        oldptr = this->release(oldIndex);
        indices_.erase(iter);
    }
    else
    {
        oldIndex = -1;
    }

    return oldptr;
}


void Foam::functionObjectList::listDir
(
    const fileName& dir,
    wordHashSet& available
)
{
    // Search specified directory for functionObject configuration files
    for (const fileName& f : fileHandler().readDir(dir))
    {
        if (f.ext().empty())
        {
            available.insert(f);
        }
    }

    // Recurse into sub-directories
    for (const fileName& d : fileHandler().readDir(dir, fileName::DIRECTORY))
    {
        listDir(dir/d, available);
    }
}


void Foam::functionObjectList::list()
{
    wordHashSet available;

    for (const fileName& d : findEtcDirs(functionObjectDictPath))
    {
        listDir(d, available);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << available.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectList::findDict(const word& funcName)
{
    // First check for functionObject dictionary file in globalCase system/

    fileName dictFile = stringOps::expand("<system>")/funcName;

    if (isFile(dictFile))
    {
        return dictFile;
    }

    for (const fileName& d : findEtcDirs(functionObjectDictPath))
    {
        dictFile = search(funcName, d);
        if (!dictFile.empty())
        {
            return dictFile;
        }
    }

    return fileName::null;
}


bool Foam::functionObjectList::readFunctionObject
(
    const string& funcNameArgs,
    dictionary& functionsDict,
    HashSet<wordRe>& requiredFields,
    const word& region
)
{
    // Parse the optional functionObject arguments:
    //     'Q(U)' -> funcName = Q; args = (U); field = U
    //
    // Supports named arguments:
    //     'patchAverage(patch=inlet, p)' -> funcName = patchAverage;
    //         args = (patch=inlet, p); field = p

    word funcName;
    wordRes args;
    List<Tuple2<word, string>> namedArgs;

    {
        const auto argsBeg = funcNameArgs.find('(');
        if (argsBeg == std::string::npos)
        {
            // Function name only, no args
            funcName = word::validate(funcNameArgs);
        }
        else
        {
            // Leading function name
            funcName = word::validate(funcNameArgs.substr(0, argsBeg));

            const auto argsEnd = funcNameArgs.rfind(')');

            stringOps::splitFunctionArgs
            (
                funcNameArgs.substr
                (
                    (argsBeg + 1),
                    (
                        (argsEnd != std::string::npos && argsBeg < argsEnd)
                      ? (argsEnd - argsBeg - 1)
                      : std::string::npos
                    )
                ),
                args,
                namedArgs
            );
        }
    }


    // Search for the functionObject dictionary
    fileName path = functionObjectList::findDict(funcName);

    if (path.empty())
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = *fileStreamPtr;

    dictionary funcsDict(fileStream);
    dictionary* funcDictPtr = funcsDict.findDict(funcName);
    dictionary& funcDict = (funcDictPtr ? *funcDictPtr : funcsDict);


    // Insert the 'field' and/or 'fields' entry corresponding to the optional
    // arguments or read the 'field' or 'fields' entry and add the required
    // fields to requiredFields
    if (args.size() == 1)
    {
        funcDict.set("field", args[0]);
        funcDict.set("fields", args);
        requiredFields.insert(args[0]);
    }
    else if (args.size() > 1)
    {
        funcDict.set("fields", args);
        requiredFields.insert(args);
    }
    else if (funcDict.found("field"))
    {
        requiredFields.insert(funcDict.get<wordRe>("field"));
    }
    else if (funcDict.found("fields"))
    {
        requiredFields.insert(funcDict.get<wordRes>("fields"));
    }

    // Insert named arguments
    for (const Tuple2<word, string>& namedArg : namedArgs)
    {
        IStringStream entryStream
        (
            namedArg.first() + ' ' + namedArg.second() + ';'
        );

        funcDict.set(entry::New(entryStream).ptr());
    }

    // Insert the region name if specified
    if (!region.empty())
    {
        funcDict.set("region", region);
    }

    // Merge this functionObject dictionary into functionsDict
    dictionary funcArgsDict;
    funcArgsDict.add(word::validate(funcNameArgs), funcDict);
    functionsDict.merge(funcArgsDict);

    return true;
}


Foam::error::handlerTypes
Foam::functionObjectList::getOrDefaultErrorHandling
(
    const word& key,
    const dictionary& dict,
    const error::handlerTypes deflt
) const
{
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);

    if (eptr)
    {
        if (eptr->isDict())
        {
            Warning
                << "The sub-dictionary '" << key
                << "' masks error handling for functions" << endl;
        }
        else
        {
            const word enumName(eptr->get<word>());

            if (!error::handlerNames.found(enumName))
            {
                // Failed the name lookup
                FatalIOErrorInFunction(dict)
                    << enumName << " is not in enumeration: "
                    << error::handlerNames << nl
                    << exit(FatalIOError);
            }

            return error::handlerNames.get(enumName);
        }
    }

    return deflt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const bool execution
)
:
    functionObjectList(runTime, runTime.controlDict(), execution)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const dictionary& parentDict,
    const bool execution
)
:
    PtrList<functionObject>(),
    errorHandling_(),
    digests_(),
    indices_(),
    warnings_(0),
    time_(runTime),
    parentDict_(parentDict),
    propsDictPtr_(nullptr),
    objectsRegistryPtr_(nullptr),
    execution_(execution),
    updated_(false)
{}


Foam::autoPtr<Foam::functionObjectList> Foam::functionObjectList::New
(
    const argList& args,
    const Time& runTime,
    dictionary& controlDict,
    HashSet<wordRe>& requiredFields
)
{
    // Merge any functions from the provided controlDict
    controlDict.add
    (
        dictionaryEntry("functions", controlDict, dictionary::null),
        true
    );

    dictionary& functionsDict = controlDict.subDict("functions");

    const word regionName = args.getOrDefault<word>("region", "");

    bool modifiedControlDict = false;

    if (args.found("dict"))
    {
        modifiedControlDict = true;

        controlDict.merge
        (
            IOdictionary
            (
                IOobject
                (
                    args["dict"],
                    runTime,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            )
        );
    }

    if (args.found("func"))
    {
        modifiedControlDict = true;

        readFunctionObject
        (
            args["func"],
            functionsDict,
            requiredFields,
            regionName
        );
    }

    if (args.found("funcs"))
    {
        modifiedControlDict = true;

        for (const word& funcName : args.getList<word>("funcs"))
        {
            readFunctionObject
            (
                funcName,
                functionsDict,
                requiredFields,
                regionName
            );
        }
    }


    autoPtr<functionObjectList> functionsPtr;

    if (modifiedControlDict)
    {
        functionsPtr.reset(new functionObjectList(runTime, controlDict));
    }
    else
    {
        functionsPtr.reset(new functionObjectList(runTime));
    }

    functionsPtr->start();

    return functionsPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::functionObjectList::triggerIndex() const
{
    return propsDict().getTrigger();
}


void Foam::functionObjectList::resetPropertiesDict()
{
    // Reset (re-read) the properties dictionary
    propsDictPtr_.reset(nullptr);
    createPropertiesDict();
}


Foam::functionObjects::properties& Foam::functionObjectList::propsDict()
{
    if (!propsDictPtr_)
    {
        createPropertiesDict();
    }

    return *propsDictPtr_;
}


const Foam::functionObjects::properties&
Foam::functionObjectList::propsDict() const
{
    if (!propsDictPtr_)
    {
        createPropertiesDict();
    }

    return *propsDictPtr_;
}


Foam::objectRegistry& Foam::functionObjectList::storedObjects()
{
    if (!objectsRegistryPtr_)
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
}


const Foam::objectRegistry& Foam::functionObjectList::storedObjects() const
{
    if (!objectsRegistryPtr_)
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
}


void Foam::functionObjectList::clear()
{
    PtrList<functionObject>::clear();
    errorHandling_.clear();
    digests_.clear();
    indices_.clear();
    warnings_.clear();
    updated_ = false;
}


Foam::label Foam::functionObjectList::findObjectID(const word& objName) const
{
    label id = 0;

    for (const functionObject& funcObj : functions())
    {
        if (funcObj.name() == objName)
        {
            return id;
        }

        ++id;
    }

    return -1;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    // For safety, also force a read() when execution is resumed
    updated_ = execution_ = false;
}


bool Foam::functionObjectList::status() const
{
    return execution_;
}


bool Foam::functionObjectList::start()
{
    return read();
}


bool Foam::functionObjectList::execute(bool writeProperties)
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        auto errIter = errorHandling_.cbegin();

        for (functionObject& funcObj : functions())
        {
            const auto errorHandling = *errIter;
            ++errIter;

            const word& objName = funcObj.name();

            if
            (
                errorHandling == error::handlerTypes::WARN
             || errorHandling == error::handlerTypes::IGNORE
            )
            {
                // Throw FatalError, FatalIOError as exceptions

                const bool oldThrowingError = FatalError.throwing(true);
                const bool oldThrowingIOerr = FatalIOError.throwing(true);

                bool hadError = false;

                // execute()
                try
                {
                    addProfiling
                    (
                        fo,
                        "functionObject::", objName, "::execute"
                    );

                    ok = funcObj.execute() && ok;
                }
                catch (const Foam::error& err)
                {
                    // Treat IOerror and error identically
                    unsigned nWarnings;
                    hadError = true;

                    if
                    (
                        (errorHandling == error::handlerTypes::WARN)
                     && (nWarnings = ++warnings_(objName)) <= maxWarnings
                    )
                    {
                        // Trickery to get original message
                        err.write(Warning, false);
                        Info<< nl
                            << "--> execute() function object '"
                            << objName << "'";

                        if (nWarnings == maxWarnings)
                        {
                            Info<< nl << "... silencing further warnings";
                        }

                        Info<< nl << endl;
                    }
                }

                if (hadError)
                {
                    // Restore previous state
                    FatalError.throwing(oldThrowingError);
                    FatalIOError.throwing(oldThrowingIOerr);
                    continue;
                }

                // write()
                try
                {
                    addProfiling
                    (
                        fo,
                        "functionObject::", objName, ":write"
                    );

                    ok = funcObj.write() && ok;
                }
                catch (const Foam::error& err)
                {
                    // Treat IOerror and error identically
                    unsigned nWarnings;
                    hadError = true;

                    if
                    (
                        (errorHandling == error::handlerTypes::WARN)
                     && (nWarnings = ++warnings_(objName)) <= maxWarnings
                    )
                    {
                        // Trickery to get original message
                        err.write(Warning, false);
                        Info<< nl
                            << "--> write() function object '"
                            << objName << "'";

                        if (nWarnings == maxWarnings)
                        {
                            Info<< nl << "... silencing further warnings";
                        }

                        Info<< nl << endl;
                    }
                }

                // Restore previous state
                FatalError.throwing(oldThrowingError);
                FatalIOError.throwing(oldThrowingIOerr);

                // Reset the warning counter (if any)
                // if no errors were encountered
                if
                (
                    (errorHandling == error::handlerTypes::WARN)
                 && !hadError && !warnings_.empty()
                )
                {
                    warnings_.erase(objName);
                }
            }
            else
            {
                // No special trapping of errors

                // execute()
                {
                    addProfiling
                    (
                        fo,
                        "functionObject::", objName, "::execute"
                    );

                    ok = funcObj.execute() && ok;
                }

                // write()
                {
                    addProfiling
                    (
                        fo,
                        "functionObject::", objName, ":write"
                    );

                    ok = funcObj.write() && ok;
                }
            }
        }
    }

    // Force writing of properties dictionary after function object execution
    if (time_.writeTime() && writeProperties)
    {
        const auto oldPrecision = IOstream::precision_;
        IOstream::precision_ = 16;

        propsDictPtr_->writeObject
        (
            IOstreamOption(IOstreamOption::ASCII, time_.writeCompression()),
            true
        );

        IOstream::precision_ = oldPrecision;
    }

    return ok;
}


bool Foam::functionObjectList::execute(const label subIndex)
{
    bool ok = execution_;

    if (ok)
    {
        for (functionObject& funcObj : functions())
        {
            // Probably do not need try/catch...

            ok = funcObj.execute(subIndex) && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::execute
(
    const UList<wordRe>& functionNames,
    const label subIndex
)
{
    bool ok = execution_;

    if (ok && functionNames.size())
    {
        for (functionObject& funcObj : functions())
        {
            if (wordRes::match(functionNames, funcObj.name()))
            {
                // Probably do not need try/catch...

                ok = funcObj.execute(subIndex) && ok;
            }
        }
    }

    return ok;
}


bool Foam::functionObjectList::end()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        auto errIter = errorHandling_.cbegin();

        for (functionObject& funcObj : functions())
        {
            const auto errorHandling = *errIter;
            ++errIter;

            const word& objName = funcObj.name();

            // Ignore failure on end() - not much we can do anyhow

            // Throw FatalError, FatalIOError as exceptions
            const bool oldThrowingError = FatalError.throwing(true);
            const bool oldThrowingIOerr = FatalIOError.throwing(true);

            try
            {
                addProfiling(fo, "functionObject::", objName, "::end");
                ok = funcObj.end() && ok;
            }
            catch (const Foam::error& err)
            {
                // Treat IOerror and error identically
                // If it somehow failed, emit a warning (unless IGNORE).
                // Unlike execute(), do not suppress further warning messages
                // (we want to know about rare issues)
                // but do reset the warnings counter for the next cycle.

                if (errorHandling != error::handlerTypes::IGNORE)
                {
                    // Trickery to get original message
                    err.write(Warning, false);
                    Info<< nl
                        << "--> end() function object '"
                        << objName << "'"
                        << nl << endl;
                }
            }

            // Restore previous state
            FatalError.throwing(oldThrowingError);
            FatalIOError.throwing(oldThrowingIOerr);

            // Reset the corresponding warning counter (if any)
            if
            (
                (errorHandling == error::handlerTypes::WARN)
             && !warnings_.empty()
            )
            {
                warnings_.erase(objName);
            }
        }
    }

    return ok;
}


bool Foam::functionObjectList::adjustTimeStep()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        for (functionObject& funcObj : functions())
        {
            const word& objName = funcObj.name();

            // Probably do not need try/catch...

            addProfiling
            (
                fo,
                "functionObject::", objName, "::adjustTimeStep"
            );

            ok = funcObj.adjustTimeStep() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::read()
{
    if (!propsDictPtr_)
    {
        createPropertiesDict();
    }

    updated_ = execution_;

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr =
        parentDict_.findEntry("functions", keyType::LITERAL);

    bool ok = true;

    if (!entryPtr)
    {
        // No functions
        PtrList<functionObject>::clear();
        errorHandling_.clear();
        digests_.clear();
        indices_.clear();
        warnings_.clear();
    }
    else if (!entryPtr->isDict())
    {
        // Bad entry type
        ok = false;
        FatalIOErrorInFunction(parentDict_)
            << "'functions' entry is not a dictionary"
            << exit(FatalIOError);
    }
    else
    {
        const dictionary& functionsDict = entryPtr->dict();

        PtrList<functionObject> newPtrs(functionsDict.size());
        List<SHA1Digest> newDigs(functionsDict.size());

        errorHandling_.resize
        (
            functionsDict.size(),
            error::handlerTypes::DEFAULT
        );

        HashTable<label> newIndices;

        addProfiling(fo, "functionObjects::read");

        // Top-level "libs" specification (optional)
        time_.libs().open
        (
            functionsDict,
            "libs",
            functionObject::dictionaryConstructorTablePtr_
        );

        // Top-level "errors" specification (optional)
        const error::handlerTypes errorHandlingFallback =
            getOrDefaultErrorHandling
            (
                "errors",
                functionsDict,
                error::handlerTypes::DEFAULT
            );

        label nFunc = 0;

        for (const entry& dEntry : functionsDict)
        {
            const word& key = dEntry.keyword();

            if (!dEntry.isDict())
            {
                // Handle or ignore some known/expected keywords

                if (key == "useNamePrefix")  // As per functionObject
                {
                    Switch sw(dEntry.stream().front());
                    if (sw.good())
                    {
                        functionObject::defaultUseNamePrefix = sw;
                    }
                    else
                    {
                        IOWarningInFunction(parentDict_)
                            << "Entry '" << key << "' is not a valid switch"
                            << endl;
                    }
                }
                else if (key != "errors" && key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry '" << key << "' is not a dictionary"
                        << endl;
                }

                continue;
            }

            const dictionary& dict = dEntry.dict();

            bool enabled = dict.getOrDefault("enabled", true);

            // Per-function "errors" specification
            const error::handlerTypes errorHandling =
                getOrDefaultErrorHandling
                (
                    "errors",
                    dict,
                    errorHandlingFallback
                );

            errorHandling_[nFunc] = errorHandling;

            newDigs[nFunc] = dict.digest();

            label oldIndex = -1;
            autoPtr<functionObject> objPtr = remove(key, oldIndex);

            const bool needsTimeControl =
                functionObjects::timeControl::entriesPresent(dict);

            if (objPtr)
            {
                // Existing functionObject:
                // Re-read if dictionary content changed and did not
                // change timeControl <-> regular

                if (enabled && newDigs[nFunc] != digests_[oldIndex])
                {
                    const bool wasTimeControl =
                        isA<functionObjects::timeControl>(*objPtr);

                    if (needsTimeControl != wasTimeControl)
                    {
                        // Changed from timeControl <-> regular

                        // Fallthrough to 'new'
                        objPtr.reset(nullptr);
                    }
                    else
                    {
                        // Normal read. Assume no errors to trap

                        addProfiling
                        (
                            fo,
                            "functionObject::", objPtr->name(), "::read"
                        );

                        enabled = objPtr->read(dict);
                    }
                }

                if (!enabled)
                {
                    // Delete disabled or an invalid(read) functionObject
                    objPtr.reset(nullptr);
                    continue;
                }
            }

            if (enabled && !objPtr)
            {
                // Throw FatalError, FatalIOError as exceptions
                const bool oldThrowingError = FatalError.throwing(true);
                const bool oldThrowingIOerr = FatalIOError.throwing(true);

                try
                {
                    // New functionObject
                    addProfiling
                    (
                        fo,
                        "functionObject::", key, "::new"
                    );
                    if (needsTimeControl)
                    {
                        objPtr.reset
                        (
                            new functionObjects::timeControl(key, time_, dict)
                        );
                    }
                    else
                    {
                        objPtr = functionObject::New(key, time_, dict);
                    }
                }
                catch (const Foam::error& err)
                {
                    objPtr.reset(nullptr);  // extra safety

                    switch (errorHandling)
                    {
                        case error::handlerTypes::IGNORE:
                            break;

                        case error::handlerTypes::STRICT:
                        {
                            exitNow(err);
                            break;
                        }

                        case error::handlerTypes::DEFAULT:
                        {
                            if (isA<Foam::IOerror>(err))
                            {
                                // Fatal for Foam::IOerror
                                exitNow(err);
                                break;
                            }

                            // Emit warning otherwise
                            [[fallthrough]];
                        }

                        case error::handlerTypes::WARN:
                        {
                            // Trickery to get original message
                            err.write(Warning, false);
                            Info<< nl
                                << "--> loading function object '"
                                << key << "'"
                                << nl << endl;
                            break;
                        }
                    }
                }

                // Restore previous state
                FatalError.throwing(oldThrowingError);
                FatalIOError.throwing(oldThrowingIOerr);

                // Require valid functionObject on all processors
                if (!returnReduceAnd(bool(objPtr)))
                {
                    objPtr.reset(nullptr);
                    ok = false;
                }
            }

            // Insert active functionObject into the list
            if (objPtr)
            {
                newPtrs.set(nFunc, std::move(objPtr));
                newIndices.insert(key, nFunc);
                ++nFunc;
            }
        }

        newPtrs.resize(nFunc);
        newDigs.resize(nFunc);
        errorHandling_.resize(nFunc);

        // Updating PtrList of functionObjects deletes any
        // existing unused functionObjects
        PtrList<functionObject>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
        warnings_.clear();
    }

    return ok;
}


bool Foam::functionObjectList::filesModified() const
{
    bool ok = false;
    if (execution_)
    {
        for (const functionObject& funcObj : functions())
        {
            bool changed = funcObj.filesModified();
            ok = ok || changed;
        }
    }
    return ok;
}


void Foam::functionObjectList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        for (functionObject& funcObj : functions())
        {
            funcObj.updateMesh(mpm);
        }
    }
}


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        for (functionObject& funcObj : functions())
        {
            funcObj.movePoints(mesh);
        }
    }
}


// ************************************************************************* //
