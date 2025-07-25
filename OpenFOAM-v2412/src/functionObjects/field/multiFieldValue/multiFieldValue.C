/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "multiFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(multiFieldValue, 0);
    addToRunTimeSelectionTable(functionObject, multiFieldValue, dictionary);
}
}
}


const Foam::Enum
<
    Foam::functionObjects::fieldValues::multiFieldValue::operationType
>
Foam::functionObjects::fieldValues::multiFieldValue::operationTypeNames_
({
    { operationType::opSum, "sum" },
    { operationType::opAdd, "add" },
    { operationType::opSubtract, "subtract" },
    { operationType::opDivide, "divide" },
    { operationType::opCmptDivide, "cmptDivide" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opAverage, "average" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "multiFieldValueImpl.C"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldValues::multiFieldValue::writeFileHeader
(
    const wordList& foNames,
    const List<wordList>& entries,
    const List<wordList>& types,
    Ostream& os
) const
{
    const word groupPrefix("Group");

    forAll(entries, i)
    {
        writeCommented(os, groupPrefix + Foam::name(i));
        os  << nl;

        forAll(entries[i], functioni)
        {
            writeCommented
            (
                os,
                "  - " + foNames[functioni] + ":" + entries[i][functioni]
            );
            os  << nl;
        }
    }

    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    forAll(entries, entryi)
    {
        writeTabbed(os, groupPrefix + Foam::name(entryi));
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::multiFieldValue::multiFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    writeFile(runTime, name, typeName, dict),
    operation_(opSubtract),
    functions_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::multiFieldValue::read
(
    const dictionary& dict
)
{
    if (!stateFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    operation_ = operationTypeNames_.get("operation", dict);

    const dictionary& functionsDict = dict.subDict("functions");
    functions_.resize(functionsDict.size());

    if (functions_.empty())
    {
        WarningInFunction
            << "No functions specified"
            << endl;
        return false;
    }

    resultFields_.resize(functions_.size());

    label functioni = 0;
    for (const entry& dEntry : functionsDict)
    {
        if (!dEntry.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Functions must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const dictionary& localDict = dEntry.dict();

        functions_.set
        (
            functioni,
            functionObject::New
            (
                IOobject::scopedName(name(), localDict.dictName()),
                time(),
                localDict
            ).ptr()
        );

        // Deactivate logging for child function objects
        //functions_[functioni].log = false;

        // Get result field names; not specified implies all
        resultFields_[functioni] =
            localDict.getOrDefault<wordList>("resultFields", wordList());

        Info<< type() << ' ' << name() << ':' << nl;
        if (resultFields_[functioni].size())
        {
            Info<< "    " << functions_[functioni].name()
                << " " << resultFields_[functioni];
        }
        else
        {
            Info<< "    " << functions_[functioni].name()
                << " - using all available entries";
        }
        Info<< nl << endl;

        ++functioni;
    }

    return true;
}


bool Foam::functionObjects::fieldValues::multiFieldValue::write()
{
    if (functions_.empty())
    {
        return false;
    }

    Log << type() << " " << name() << " write:" << endl;

    const label nFunction = functions_.size();
    wordList entries0;
    label nEntries = -1;

    wordList foNames(nFunction);
    List<wordList> entries;
    List<wordList> types;

    forAll(functions_, functioni)
    {
        auto& f = functions_[functioni];
        foNames[functioni] = f.name();

        // Note: replicating functionObjectList execute() and write()
        // - results may be written on either
        f.execute();
        f.write();

        wordList e = resultFields_[functioni];
        if (e.empty())
        {
            e = objectResultEntries(f.name());
        }

        if (functioni == 0)
        {
            entries0 = e;
            nEntries = e.size();
            entries.resize(nEntries);
            types.resize(nEntries);

            forAll(entries, entryi)
            {
                entries[entryi].resize(nFunction);
                types[entryi].resize(nFunction);
            }
        }

        if (e.size() != nEntries)
        {
            const word& f0Name = functions_[0].name();

            FatalErrorInFunction
                << "Inconsistent number of result entries" << nl
                << "    " << f0Name << " entries:" << entries0 << nl
                << "    " << f.name() << " entries:" << e
                << exit(FatalError);
        }

        forAll(e, entryi)
        {
            entries[entryi][functioni] = e[entryi];
            types[entryi][functioni] = objectResultType(f.name(), e[entryi]);

            if (types[entryi][functioni] == word::null)
            {
                FatalErrorInFunction
                    << "Unable to find function object result" << nl
                    << "    function object   : " << f.name() << nl
                    << "    result name       : " << e[entryi] << nl
                    << "    available results : "
                        << objectResultEntries(f.name())
                    << exit(FatalError);
            }
        }
    }

    if (!writtenHeader_)
    {
        writeFileHeader(foNames, entries, types, file());
        writtenHeader_ = true;
    }

    writeCurrentTime(file());

    forAll(entries, i)
    {
        const wordList& entryi = entries[i];
        const word& expectedType = types[i][0];
        const wordList& foTypes = types[i];

        forAll(foTypes, functioni)
        {
            const word& foType = foTypes[functioni];

            if (foType != expectedType)
            {
                FatalErrorInFunction
                    << "Inconsistent function result types" << nl
                    << "    " << functions_[0].name()
                    << " result type:" << expectedType << nl
                    << "    " << functions_[functioni].name()
                    << " result type:" << foType
                    << exit(FatalError);
            }
        }

        const bool ok
        (
            applyOperation<scalar>(expectedType, foNames, entryi)
         || applyOperation<vector>(expectedType, foNames, entryi)
         || applyOperation<sphericalTensor>(expectedType, foNames, entryi)
         || applyOperation<symmTensor>(expectedType, foNames, entryi)
         || applyOperation<tensor>(expectedType, foNames, entryi)
        );

        if (!ok)
        {
            Log << "Operation not applied between functions:" << nl
                << flatOutput(foNames, FlatOutput::BareComma{}) << nl
                << "with result names:" << nl
                << flatOutput(entryi, FlatOutput::BareComma{})
                << endl;
        }
    }

    Log << (nEntries == 0 ? "    none" : "") << endl;

    file()<< endl;

    return true;
}


bool Foam::functionObjects::fieldValues::multiFieldValue::execute()
{
    return true;
}


// ************************************************************************* //
