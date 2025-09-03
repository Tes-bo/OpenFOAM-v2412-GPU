/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    foamDictionary

Description
    Interrogate and manipulate dictionaries.

Usage
    \b foamDictionary [OPTION] dictionary

      - \par -entry \<name\>
        Selects an entry

      - \par -keywords
        Prints the keywords (of the selected entry or of the top level if
        no entry was selected

      - \par -add \<value\>
        Adds the entry (should not exist yet)

      - \par -set \<value\>
        Adds or replaces the entry

      - \par -remove
        Remove the selected entry

      - \par -diff \<dictionary\>
        Write differences with respect to the specified dictionary
        (or sub entry if -entry specified)

      - \par -diff-etc \<dictionary\>
        Write differences with respect to the specified dictionary
        (or sub entry if -entry specified)

      - \par -expand
        Read the specified dictionary file, expand the macros etc. and write
        the resulting dictionary to standard output.

      - \par -includes
        List the \c \#include and \c \#sinclude files to standard output

      - \par -disableFunctionEntries
        Do not expand macros or directives (\#include etc)

      - \par -precision int
        Set default write precision for IOstreams

    Example usage:
      - Change simulation to run for one timestep only:
        \verbatim
          foamDictionary system/controlDict -entry stopAt -set writeNow
        \endverbatim

      - Change solver:
        \verbatim
           foamDictionary system/fvSolution -entry solvers/p/solver -set PCG
        \endverbatim

      - Print bc type:
        \verbatim
           foamDictionary 0/U -entry boundaryField/movingWall/type
        \endverbatim

      - Change bc parameter:
        \verbatim
           foamDictionary 0/U -entry boundaryField/movingWall/value \
             -set "uniform (2 0 0)"
        \endverbatim

      - Change whole bc type:
        \verbatim
          foamDictionary 0/U -entry boundaryField/movingWall \
            -set "{type uniformFixedValue; uniformValue (2 0 0);}"
        \endverbatim

      - Write the differences with respect to a template dictionary:
        \verbatim
          foamDictionary 0/U -diff-etc templates/closedVolume/0/U
        \endverbatim

      - Write the differences in boundaryField with respect to a
        template dictionary:
        \verbatim
          foamDictionary 0/U -diff-etc templates/closedVolume/0/U \
            -entry boundaryField
        \endverbatim

      - Change patch type:
        \verbatim
          foamDictionary constant/polyMesh/boundary \
            -entry entry0/fixedWalls/type -set patch
        \endverbatim
        This uses special parsing of Lists which stores these in the
        dictionary with keyword 'entryDDD' where DDD is the position
        in the dictionary (after ignoring the FoamFile entry).

    Notes:
        - the use of '.' as the scoping symbol might conflict with
        e.g. file extensions ('.' is not really considered
        to be a valid word character). Instead use the '/' as a scoping
        character e.g.
          foamDictionary system/snappyHexMeshDict \
            -entry /geometry/motorBike.obj -remove

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "Time.H"
#include "Fstream.H"
#include "etcFiles.H"
#include "includeEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Convert very old ':' scope syntax to less old '.' scope syntax,
//  but leave anything with '/' delimiters untouched
bool upgradeScope(word& entryName)
{
    if (entryName.contains(':') && !entryName.contains('/'))
    {
        InfoErr
            << "Warning: upgrading very old ':' scope syntax: \""
            << entryName << '"' << endl;

        // Make copy - cannot use stringOps::split
        const wordList cmpts(fileName(entryName).components(':'));

        entryName.clear();

        bool addSep = false;

        for (const word& cmpt : cmpts)
        {
            if (addSep) entryName += '.';
            if (!cmpt.empty())
            {
                addSep = true;
                entryName += cmpt;
            }
        }

        return true;
    }

    // Nothing changed
    return false;
}


//- Split into dictionary name and the entry name
class dictAndKeyword
{
    word dict_;
    word key_;

public:

    dictAndKeyword(const word& scopedName)
    {
        auto i = scopedName.rfind('/');
        if (i == std::string::npos)
        {
            i = scopedName.rfind('.');
        }

        if (i != std::string::npos)
        {
            dict_ = scopedName.substr(0, i);
            key_  = scopedName.substr(i+1);
        }
        else
        {
            key_  = scopedName;
        }
    }

    const word& dict() const noexcept { return dict_; }

    const word& key() const noexcept { return key_; }
};


const dictionary& lookupScopedDict
(
    const dictionary& dict,
    const word& subDictName
)
{
    if (subDictName.empty())
    {
        return dict;
    }

    const entry* eptr = dict.findScoped(subDictName, keyType::LITERAL);

    if (!eptr || !eptr->isDict())
    {
        FatalIOErrorInFunction(dict)
            << "'" << subDictName << "' not found in dictionary "
            << dict.name() << " or is not a dictionary" << nl
            << "Known entries are " << dict.keys()
            << exit(FatalIOError);
    }

    return eptr->dict();
}


void removeDict(dictionary& dict, const dictionary& dictToRemove)
{
    for (const entry& refEntry : dictToRemove)
    {
        auto finder = dict.search(refEntry.keyword(), keyType::LITERAL);

        bool purge = false;

        if (finder.isDict())
        {
            if (refEntry.isDict())
            {
                removeDict(finder.dict(), refEntry.dict());

                // Purge if dictionary is empty
                purge = finder.dict().empty();
            }
        }
        else if (finder.good() && !refEntry.isDict())
        {
            // Purge if entries match
            purge = (finder.ref() == refEntry);
        }

        if (purge)
        {
            dict.remove(refEntry.keyword());
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Interrogate and manipulate dictionaries"
    );

    argList::noBanner();
    argList::noJobInfo();
    argList::addArgument("dict", "The dictionary file to process");
    argList::addBoolOption("keywords", "List keywords");
    argList::addOption("entry", "name", "Report/select the named entry");
    argList::addBoolOption
    (
        "value",
        "Print entry value"
    );
    argList::addOption
    (
        "set",
        "value",
        "Set entry value or add new entry"
    );
    argList::addOption
    (
        "add",
        "value",
        "Add a new entry"
    );
    argList::addBoolOption
    (
        "remove",
        "Remove the entry"
    );
    argList::addOption
    (
        "diff",
        "dict",
        "Write differences with respect to the specified dictionary"
    );
    argList::addOption
    (
        "diff-etc",
        "dict",
        "As per -diff, but locate the file as per foamEtcFile"
    );
    argList::addOptionCompat("diff-etc", {"diffEtc", 1712});
    argList::addOption
    (
        "precision",
        "int",
        "Set default write precision for IOstreams"
    );

    argList::addBoolOption
    (
        "includes",
        "List the #include/#sinclude files to standard output"
    );
    argList::addBoolOption
    (
        "expand",
        "Read the specified dictionary file, expand the macros etc. and write "
        "the resulting dictionary to standard output"
    );
    argList::addBoolOption
    (
        "disableFunctionEntries",
        "Disable expansion of dictionary directives - #include, #codeStream etc"
    );
    profiling::disable(); // Disable profiling (and its output)

    argList args(argc, argv);

    const bool listIncludes = args.found("includes");

    if (listIncludes)
    {
        Foam::functionEntries::includeEntry::log = true;
    }

    const bool disableEntries = args.found("disableFunctionEntries");
    if (disableEntries)
    {
        // Report on stderr (once) to avoid polluting the output
        InfoErr<< "Not expanding variables or dictionary directives" << endl;
        entry::disableFunctionEntries = true;
    }

    // Set the default output precision
    {
        const unsigned prec = args.getOrDefault<unsigned>("precision", 0u);
        if (prec)
        {
            // InfoErr<< "Output write precision set to " << prec << endl;
            IOstream::defaultPrecision(prec);
            Sout.precision(prec);
        }
    }

    const auto dictFileName = args.get<fileName>(1);

    auto dictFile = autoPtr<IFstream>::New(dictFileName);
    if (!dictFile().good())
    {
        FatalErrorInFunction
            << "Cannot open file " << dictFileName
            << exit(FatalError, 1);
    }


    bool changed = false;

    // Read but preserve headers
    dictionary dict(dictFile(), true);

    if (listIncludes)
    {
        return 0;
    }
    else if (args.found("expand"))
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictFileName << "\n//\n";
        dict.write(Info, false);
        IOobject::writeDivider(Info);

        return 0;
    }


    // Has "diff" or "diff-etc"
    bool optDiff = false;

    // Reference dictionary for -diff / -diff-etc
    dictionary diffDict;
    {
        fileName diffFileName;
        if (args.readIfPresent("diff", diffFileName))
        {
            IFstream diffFile(diffFileName);
            if (!diffFile.good())
            {
                FatalErrorInFunction
                    << "Cannot open file " << diffFileName
                    << exit(FatalError, 1);
            }

            // Read but preserve headers
            diffDict.read(diffFile, true);
            optDiff = true;
        }
        else if (args.readIfPresent("diff-etc", diffFileName))
        {
            fileName foundName = findEtcFile(diffFileName);
            if (foundName.empty())
            {
                FatalErrorInFunction
                    << "Cannot find etcFile " << diffFileName
                    << exit(FatalError, 1);
            }

            IFstream diffFile(foundName);
            if (!diffFile.good())
            {
                FatalErrorInFunction
                    << "Cannot open file " << foundName
                    << exit(FatalError, 1);
            }

            // Read but preserve headers
            diffDict.read(diffFile, true);
            optDiff = true;
        }
    }

    word scopedName;  // Actually fileName, since it can contain '/' scoping
    if (args.readIfPresent("entry", scopedName))
    {
        upgradeScope(scopedName);

        string newValue;
        if
        (
            args.readIfPresent("set", newValue)
         || args.readIfPresent("add", newValue)
        )
        {
            const bool overwrite = args.found("set");

            // Dictionary name and keyword
            const dictAndKeyword dAk(scopedName);

            // The context for the action
            const dictionary& d(lookupScopedDict(dict, dAk.dict()));

            // Create a new entry
            IStringStream str(string(dAk.key()) + ' ' + newValue + ';');
            entry* ePtr(entry::New(str).ptr());

            if (overwrite)
            {
                const_cast<dictionary&>(d).set(ePtr);
            }
            else
            {
                const_cast<dictionary&>(d).add(ePtr, false);
            }
            changed = true;

            const auto finder = dict.csearchScoped(scopedName, keyType::REGEX);

            // Print the changed entry to stderr
            if (finder.good())
            {
                InfoErr<< finder.ref();
            }
        }
        else if (args.found("remove"))
        {
            // Dictionary name and keyword
            const dictAndKeyword dAk(scopedName);

            // The context for the action
            const dictionary& d(lookupScopedDict(dict, dAk.dict()));

            const_cast<dictionary&>(d).remove(dAk.key());
            changed = true;
        }
        else
        {
            // Optionally remove a second dictionary
            if (optDiff)
            {
                // Dictionary name and keyword
                const dictAndKeyword dAk(scopedName);

                const dictionary& d1(lookupScopedDict(dict, dAk.dict()));
                const dictionary& d2(lookupScopedDict(diffDict, dAk.dict()));

                const entry* e1Ptr = d1.findEntry(dAk.key(), keyType::REGEX);
                const entry* e2Ptr = d2.findEntry(dAk.key(), keyType::REGEX);

                if (e1Ptr && e2Ptr)
                {
                    if (*e1Ptr == *e2Ptr)
                    {
                        const_cast<dictionary&>(d1).remove(dAk.key());
                    }
                    else if (e1Ptr->isDict() && e2Ptr->isDict())
                    {
                        removeDict
                        (
                            const_cast<dictionary&>(e1Ptr->dict()),
                            e2Ptr->dict()
                        );
                    }
                }
            }

            const auto finder = dict.csearchScoped(scopedName, keyType::REGEX);

            if (!finder.good())
            {
                FatalIOErrorInFunction(dictFile())
                    << "Cannot find entry " << scopedName
                    << exit(FatalIOError, 2);
            }
            else if (args.found("keywords"))
            {
                // Report keywords to stdout
                for (const entry& e : finder.dict())
                {
                    Info<< e.keyword() << endl;
                }
            }
            else if (args.found("value"))
            {
                // Report value to stdout
                if (finder.isDict())
                {
                    Info<< finder.dict();
                }
                else if (finder.isStream())
                {
                    bool separator = false;

                    for (const token& tok : finder.stream())
                    {
                        if (separator)
                        {
                            Info<< token::SPACE;
                        }
                        separator = true;
                        Info<< tok;
                    }
                    Info<< endl;
                }
            }
            else
            {
                // Report entry to stdout
                Info<< finder.ref();
            }
        }
    }
    else if (args.found("keywords"))
    {
        // Report keywords to stdout
        for (const entry& e : dict)
        {
            Info<< e.keyword() << endl;
        }
    }
    else if (optDiff)
    {
        // Report difference to stdout
        removeDict(dict, diffDict);
        dict.write(Info, false);
    }
    else
    {
        // Report dictionary to stdout
        dict.write(Info, false);
    }

    if (changed)
    {
        dictFile.clear();
        OFstream os(dictFileName);
        IOobject::writeBanner(os);
        IOobject::writeDivider(os);
        dict.write(os, false);
        IOobject::writeEndDivider(os);
    }

    return 0;
}


// ************************************************************************* //
