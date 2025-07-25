/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    setExprBoundaryFields

Group
    grpPreProcessingUtilities

Description
    Set boundary values using an expression

Note
    Based on funkySetBoundaryFields from
    Bernhard Gschaider <bgschaid@hfd-research.com>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "patchExprDriver.H"
#include "timeSelector.H"
#include "readFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Normally without functionObjects, with -withFunctionObjects to enable
    argList::noFunctionObjects(true);

    // No -constant, no special treatment for 0/
    timeSelector::addOptions(false);

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of the controlDict setting"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for setExprBoundaryFieldsDict"
    );
    argList::addBoolOption
    (
        "cache-fields",
        "Cache fields between calls",
        true // Advanced
    );
    argList::addOption
    (
        "load-fields",
        "wordList",
        "Specify field or fields to preload. Eg, 'T' or '(p T U)'",
        true // Advanced
    );
    argList::addBoolOption
    (
        "backup",
        "Preserve sub-entry as .backup",
        true // Advanced
    );
    argList::addDryRunOption
    (
        "Evaluate but do not write"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"

    const bool dryrun      = args.dryRun();
    const bool backup      = args.found("backup");
    const bool cacheFields = args.found("cache-fields");

    if (cacheFields)
    {
        Warning
            << "The current cache-fields behaviour (caching disk reads) "
            << "may lead to unexpected behaviour as previous modifications "
            << "will not be visible."
            << endl;
    }

    const word dictName("setExprBoundaryFieldsDict");

    #include "createTime.H"

    instantList times = timeSelector::select0(runTime, args);

    if (times.empty())
    {
        FatalErrorInFunction
            << "No times selected." << nl
            << exit(FatalError);
    }

    #include "createNamedMesh.H"

    #include "setSystemMeshDictionaryIO.H"
    IOdictionary setExprDict(dictIO);

    IOstreamOption streamOpt(runTime.writeFormat());
    if (args.found("ascii"))
    {
        streamOpt.format(IOstreamOption::ASCII);
    }

    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "\nTime = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // preload fields specified on command-line
        if (timei == 0)
        {
            wordList preloadFields;
            args.readListIfPresent("load-fields", preloadFields);
            readFieldsHandler(mesh).execute(preloadFields);
        }
        // preload fields specified in dictionary
        {
            wordList preloadFields;
            setExprDict.readIfPresent("readFields", preloadFields);
            readFieldsHandler(mesh).execute(preloadFields);
        }

        for (const entry& dEntry : setExprDict)
        {
            if (!dEntry.isDict())
            {
                if (dEntry.keyword() != "readFields")
                {
                    Info<< "Ignoring non-dictionary entry "
                        << dEntry.keyword() << nl;
                }
                continue;
            }

            const dictionary& dict = dEntry.dict();

            const word fieldName(dict.get<word>("field"));

            List<dictionary> exprDicts;
            dict.readEntry("expressions", exprDicts);

            if (exprDicts.empty())
            {
                Info<< "No expressions for " << fieldName << nl;
                continue;
            }


            // Read dictionary
            // Note: disable class type checking so we can load field
            const word oldTypeName = IOdictionary::typeName;
            const_cast<word&>(IOdictionary::typeName) = word::null;

            IOobject fieldHeader
            (
                fieldName,
                mesh.thisDb().time().timeName(),
                mesh.thisDb(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            const bool headOk = fieldHeader.typeHeaderOk<IOdictionary>(false);

            if (!headOk)
            {
                // Restore type
                const_cast<word&>(IOdictionary::typeName) = oldTypeName;

                WarningInFunction
                    << "Requested field to change " << fieldName
                    << " does not exist in " << fieldHeader.path() << endl;
                continue;
            }

            IOdictionary fieldDict(fieldHeader);

            // Restore type
            const_cast<word&>(IOdictionary::typeName) = oldTypeName;

            // Fake type back to what was in field
            const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();

            Info<< "Processing field " << fieldName << nl;

            dictionary& boundaryFieldDict = fieldDict.subDict("boundaryField");

            for (const dictionary& currDict : exprDicts)
            {
                const word patchName(currDict.get<word>("patch"));
                const word targetName(currDict.get<word>("target"));

                dictionary& patchDict = boundaryFieldDict.subDict(patchName);

                expressions::exprString valueExpr_("expression", currDict);

                Info<< "Set boundaryField/" << patchName << '/'
                    << targetName << nl
                    << "with expression" << nl
                    << "<<<<" << nl
                    << valueExpr_.c_str() << nl
                    << ">>>>" << nl;

                expressions::patchExprDriver driver(currDict, mesh);

                // Search files only
                driver.setSearchBehaviour
                (
                    expressions::exprDriver::SEARCH_FILES,
                    cacheFields
                );

                driver.clearVariables();
                driver.parse(valueExpr_);

                // Serializing via Field::writeEntry etc
                OStringStream serialize;
                driver.result().writeEntry("", serialize);

                if (backup && !dryrun)
                {
                    patchDict.changeKeyword
                    (
                        targetName,
                        word(targetName + ".backup"),
                        true  // Overwrite
                    );
                }

                patchDict.set(targetName, serialize.str().c_str());

                if (dryrun)
                {
                    Info<< "Evaluated:" << nl
                        << "<<<<" << nl
                        << serialize.str().c_str()  // (already includes nl)
                        << ">>>>" << nl;
                }
            }

            if (!dryrun)
            {
                Info<< "Write " << fieldDict.filePath() << nl;
                fieldDict.regIOobject::writeObject(streamOpt, true);
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
