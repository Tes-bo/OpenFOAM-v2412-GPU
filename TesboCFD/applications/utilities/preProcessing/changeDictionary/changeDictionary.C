/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Application
    changeDictionary

Group
    grpPreProcessingUtilities

Description
    Utility to change dictionary entries, e.g. can be used to change the patch
    type in the field and polyMesh/boundary files.

    Reads dictionaries (fields) and entries to change from a dictionary.
    E.g. to make the \em movingWall a \em fixedValue for \em p but all other
    \em Walls a zeroGradient boundary condition, the
    \c system/changeDictionaryDict would contain the following:
    \verbatim
    p                           // field to change
    {
        boundaryField
        {
            ".*Wall"            // entry to change
            {
                type            zeroGradient;
            }
            movingWall          // entry to change
            {
                type            fixedValue;
                value           uniform 123.45;
            }
        }
    }
    \endverbatim
    Replacement entries starting with '~' will remove the entry.

Usage
    \b changeDictionary [OPTION]

    Options:
      - \par -subDict
        Specify the subDict name of the replacements dictionary.

      - \par -literalRE
        Do not interpret regular expressions or patchGroups; treat them as any
        other keyword.

      - \par -enableFunctionEntries
        Enable function entries (default: disabled)

      - \par -disablePatchGroups
        Disable the default checking for keys being patchGroups

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Extract groupPatch info from boundary file info
HashTable<wordList> extractPatchGroups(const dictionary& boundaryDict)
{
    HashTable<wordList> groupToPatch;

    for (const entry& dEntry : boundaryDict)
    {
        if (!dEntry.isDict())
        {
            continue;
        }

        const word& patchName = dEntry.keyword();
        const dictionary& patchDict = dEntry.dict();

        wordList groupNames;
        patchDict.readIfPresent("inGroups", groupNames);

        for (const word& groupName : groupNames)
        {
            auto groupIter = groupToPatch.find(groupName);
            if (groupIter.good())
            {
                (*groupIter).append(patchName);
            }
            else
            {
                groupToPatch.insert(groupName, wordList(one{}, patchName));
            }
        }
    }

    return groupToPatch;
}


bool merge
(
    const bool addNonExisting,
    dictionary&,
    const dictionary&,
    const bool,
    const HashTable<wordList>&
);


// Add thisEntry to dictionary thisDict.
bool addEntry
(
    dictionary& thisDict,
    entry& thisEntry,
    const entry& mergeEntry,
    const bool literalRE,
    const HashTable<wordList>& shortcuts
)
{
    bool changed = false;

    // Recursively merge sub-dictionaries
    // TODO: merge without copying
    if (thisEntry.isDict() && mergeEntry.isDict())
    {
        if
        (
            merge
            (
                true,
                const_cast<dictionary&>(thisEntry.dict()),
                mergeEntry.dict(),
                literalRE,
                shortcuts
            )
        )
        {
            changed = true;
        }
    }
    else
    {
        // Should use in-place modification instead of adding
        thisDict.add(mergeEntry.clone(thisDict).ptr(), true);
        changed = true;
    }

    return changed;
}


// List of indices into thisKeys
labelList findMatches
(
    const HashTable<wordList>& shortcuts,
    const wordList& shortcutNames,
    const wordList& thisKeys,
    const wordRe& key
)
{
    labelList matches;

    if (key.isPattern())
    {
        // Wildcard match
        matches = wordRes::matching(key, thisKeys);
    }
    else if (shortcuts.size())
    {
        // See if patchGroups expand to valid thisKeys
        labelList indices = wordRes::matching(key, shortcutNames);

        for (const label idx : indices)
        {
            const word& name = shortcutNames[idx];
            const wordList& keys = shortcuts[name];
            for (const word& k : keys)
            {
                const label index = thisKeys.find(k);
                if (index != -1)
                {
                    matches.append(index);
                }
            }
        }
    }
    return matches;
}


// Dictionary merging/editing.
// literalRE:
// - true: behave like dictionary::merge, i.e. add regexps just like
//   any other key.
// - false : interpret wildcard as a rule for items to be matched.
bool merge
(
    const bool addNonExisting,
    dictionary& thisDict,
    const dictionary& mergeDict,
    const bool literalRE,
    const HashTable<wordList>& shortcuts
)
{
    const wordList shortcutNames(shortcuts.toc());

    bool changed = false;

    // Save current (non-wildcard) keys before adding items.
    wordHashSet thisKeysSet;
    {
        for (const word& k : thisDict.keys(false))
        {
            thisKeysSet.insert(k);
        }
    }

    // Pass 1. All literal matches

    for (const entry& mergeEntry : mergeDict)
    {
        const keyType& key = mergeEntry.keyword();

        if (key[0] == '~')
        {
            const word eraseKey = key.substr(1);
            if (thisDict.remove(eraseKey))
            {
                // Mark thisDict entry as having been match for wildcard
                // handling later on.
                thisKeysSet.erase(eraseKey);
            }
            changed = true;
        }
        else if (literalRE || !(key.isPattern() || shortcuts.found(key)))
        {
            entry* eptr = thisDict.findEntry(key, keyType::LITERAL);

            if (eptr)
            {
                // Mark thisDict entry as having been match for wildcard
                // handling later on.
                thisKeysSet.erase(eptr->keyword());

                if
                (
                    addEntry
                    (
                        thisDict,
                       *eptr,
                        mergeEntry,
                        literalRE,
                        shortcuts
                    )
                )
                {
                    changed = true;
                }
            }
            else
            {
                if (addNonExisting)
                {
                    // Not found - just add
                    thisDict.add(mergeEntry.clone(thisDict).ptr());
                    changed = true;
                }
                else
                {
                    IOWarningInFunction(mergeDict)
                        << "Ignoring non-existing entry " << key
                        << endl;
                }
            }
        }
    }


    // Pass 2. Wildcard or shortcut matches (if any) on any non-match keys.

    if (!literalRE && thisKeysSet.size())
    {
        // Pick up remaining dictionary entries
        wordList thisKeys(thisKeysSet.toc());

        for (const entry& mergeEntry : mergeDict)
        {
            const keyType& key = mergeEntry.keyword();

            if (key[0] == '~')
            {
                const word eraseKey = key.substr(1);

                // List of indices into thisKeys
                labelList matches
                (
                    findMatches
                    (
                        shortcuts,
                        shortcutNames,
                        thisKeys,
                        eraseKey
                    )
                );

                // Remove all matches
                for (const label matchi : matches)
                {
                    const word& k = thisKeys[matchi];
                    thisKeysSet.erase(k);
                }
                changed = true;
            }
            else
            {
                // List of indices into thisKeys
                labelList matches
                (
                    findMatches
                    (
                        shortcuts,
                        shortcutNames,
                        thisKeys,
                        key
                    )
                );

                // Add all matches
                for (const label matchi : matches)
                {
                    const word& k = thisKeys[matchi];

                    entry* eptr = thisDict.findEntry(k, keyType::LITERAL);

                    if
                    (
                        addEntry
                        (
                            thisDict,
                           *eptr,
                            mergeEntry,
                            literalRE,
                            HashTable<wordList>(0)    // no shortcuts
                                                      // at deeper levels
                        )
                    )
                    {
                        changed = true;
                    }
                }
            }
        }
    }

    return changed;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Utility to change dictionary entries"
        " (such as the patch type for fields and polyMesh/boundary files)."
    );

    argList::addOption("dict", "file", "Alternative changeDictionaryDict");

    argList::addOption
    (
        "subDict",
        "name",
        "Specify the subDict name of the replacements dictionary"
    );
    argList::addOption
    (
        "instance",
        "name",
        "Override instance setting (default is the time name)"
    );

    // Add explicit time option
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "literalRE",
        "Treat regular expressions literally (i.e., as a keyword)"
    );
    argList::addBoolOption
    (
        "enableFunctionEntries",
        "Enable expansion of dictionary directives - #include, #codeStream etc"
    );
    argList::addBoolOption
    (
        "disablePatchGroups",
        "Disable matching keys to patch groups"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    // Optionally override controlDict time with -time options
    instantList times = timeSelector::selectIfPresent(runTime, args);
    if (times.size() < 1)
    {
        FatalErrorInFunction
            << "No times selected." << exit(FatalError);
    }
    forAll(times, timei)
    {
        word instance;
        if (args.readIfPresent("instance", instance))
        {
            if (times.size() > 1)
            {
                FatalErrorInFunction
                    << "Multiple times selected with 'instance' option"
                    << exit(FatalError);
            }
        }
        else
        {
            runTime.setTime(times[timei], timei);
            instance = runTime.timeName();
        }

        #include "createNamedMesh.H"

        const bool literalRE = args.found("literalRE");
        if (literalRE)
        {
            Info<< "Not interpreting any regular expressions (RE)"
                << " in the changeDictionaryDict." << endl
                << "Instead they are handled as any other entry, i.e. added if"
                << " not present." << endl;
        }

        const bool enableEntries = args.found("enableFunctionEntries");
        if (enableEntries)
        {
            Info<< "Allowing dictionary preprocessing (#include, #codeStream)."
                << endl;
        }

        const int oldFlag = entry::disableFunctionEntries;
        if (!enableEntries)
        {
            // By default disable dictionary expansion for fields
            entry::disableFunctionEntries = 1;
        }


        const bool disablePatchGroups = args.found("disablePatchGroups");
        if (disablePatchGroups)
        {
            Info<< "Not interpreting any keys in the changeDictionary"
                << " as patchGroups"
                << endl;
        }


        // Make sure we do not use the master-only reading since we read
        // fields (different per processor) as dictionaries.
        IOobject::fileModificationChecking = IOobject::timeStamp;


        // Get the replacement rules from a dictionary

        const word dictName("changeDictionaryDict");
        #include "setSystemMeshDictionaryIO.H"
        IOdictionary dict(dictIO);

        const dictionary* replaceDictsPtr = &dict;

        if (args.found("subDict"))
        {
            replaceDictsPtr = &dict.subDict(args["subDict"]);
        }

        const dictionary& replaceDicts = *replaceDictsPtr;

        Info<< "Read dictionary " << dict.name()
            << " with replacements for dictionaries "
            << replaceDicts.toc() << endl;



        // Always read boundary to get patch groups
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Reading polyMesh/boundary file to extract patch names"
            << endl;

        // Read PtrList of dictionary as dictionary.
        const word oldTypeName = IOPtrList<entry>::typeName;
        const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
        IOPtrList<entry> dictList
        (
            IOobject
            (
                "boundary",
                runTime.findInstance
                (
                    polyMesh::meshDir(regionName),
                    "boundary",
                    IOobject::READ_IF_PRESENT
                ),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;

        // Fake type back to what was in field
        const_cast<word&>(dictList.type()) = dictList.headerClassName();

        // Temporary convert to dictionary
        dictionary fieldDict;
        for (const entry& e : dictList)
        {
            if (e.isDict())
            {
                fieldDict.add(e.keyword(), e.dict());
            }
        }

        if (dictList.size())
        {
            Info<< "Loaded dictionary " << dictList.name()
                << " with entries " << fieldDict.toc() << endl;
        }

        // Extract any patchGroups information (= shortcut for set of
        // patches)
        HashTable<wordList> patchGroups;
        if (!disablePatchGroups)
        {
            patchGroups = extractPatchGroups(fieldDict);
            if (patchGroups.size())
            {
                Info<< "Extracted patch groups:" << endl;
                wordList groups(patchGroups.sortedToc());
                forAll(groups, i)
                {
                    Info<< "    group " << groups[i] << " with patches "
                        << patchGroups[groups[i]] << endl;
                }
            }
        }


        // Every replacement is a dictionary name and a keyword in this

        for (const entry& replaceEntry : replaceDicts)
        {
            if (!replaceEntry.isDict())
            {
                // Could also warn
                continue;
            }

            const word& fieldName = replaceEntry.keyword();
            const dictionary& replaceDict = replaceEntry.dict();

            Info<< "Replacing entries in dictionary " << fieldName << endl;

            // Handle 'boundary' specially:
            // - is PtrList of dictionaries
            // - is in polyMesh/
            if (fieldName == "boundary")
            {
                Info<< "Special handling of " << fieldName
                    << " as polyMesh/boundary file." << endl;

                // Merge the replacements in. Do not add non-existing entries.
                Info<< "Merging entries from " << replaceDict.toc() << endl;
                merge(false, fieldDict, replaceDict, literalRE, patchGroups);

                Info<< "fieldDict:" << fieldDict << endl;

                // Convert back into dictList
                wordList doneKeys(dictList.size());

                label nEntries = fieldDict.size();
                nEntries = 0;

                forAll(dictList, i)
                {
                    doneKeys[i] = dictList[i].keyword();

                    const entry* ePtr = fieldDict.findEntry
                    (
                        doneKeys[i],
                        keyType::REGEX
                    );
                    // Check that it hasn't been removed from fieldDict
                    if (ePtr)
                    {
                        dictList.set(nEntries++, ePtr->clone());
                        fieldDict.remove(doneKeys[i]);
                    }
                }

                // Add remaining entries
                for (const entry& e : fieldDict)
                {
                    dictList.set(nEntries++, e.clone());
                }
                dictList.setSize(nEntries);

                Info<< "Writing modified " << fieldName << endl;
                dictList.writeObject
                (
                    IOstreamOption(runTime.writeFormat()),
                    true
                );
            }
            else
            {
                // Read dictionary
                // Note: disable class type checking so we can load field
                Info<< "Loading dictionary " << fieldName << endl;
                const word oldTypeName = localIOdictionary::typeName;
                const_cast<word&>(localIOdictionary::typeName) = word::null;

                IOobject fieldHeader
                (
                    fieldName,
                    instance,
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                if (fieldHeader.typeHeaderOk<localIOdictionary>(false))
                {
                    //IOdictionary fieldDict(fieldHeader);
                    //- dictionaries to-be-changed are either boundary
                    //  or field dictionary. Both are processor-local.
                    localIOdictionary fieldDict(fieldHeader);

                    const_cast<word&>(localIOdictionary::typeName) =
                        oldTypeName;

                    // Fake type back to what was in field
                    const_cast<word&>(fieldDict.type()) =
                        fieldDict.headerClassName();

                    Info<< "Loaded dictionary " << fieldName
                        << " with entries " << fieldDict.toc() << endl;

                    // Merge the replacements in (allow adding)
                    Info<< "Merging entries from " << replaceDict.toc() << endl;
                    merge(true, fieldDict, replaceDict, literalRE, patchGroups);

                    Info<< "Writing modified fieldDict " << fieldName << endl;
                    fieldDict.regIOobject::write();
                }
                else
                {
                    WarningInFunction
                        << "Requested field to change " << fieldName
                        << " does not exist in " << fieldHeader.path() << endl;
                }
            }

            entry::disableFunctionEntries = oldFlag;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
