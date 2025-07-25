/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    foamUpgradeCyclics

Group
    grpPreProcessingUtilities

Description
    Tool to upgrade mesh and fields for split cyclics.

Usage
    \b foamUpgradeCyclics [OPTION]

    Options:
      - \par -dry-run
        Suppress writing the updated files with split cyclics

      - \par -enableFunctionEntries
        By default all dictionary preprocessing of fields is disabled

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "entry.H"
#include "IOPtrList.H"
#include "cyclicPolyPatch.H"
#include "dictionaryEntry.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// Read boundary file without reading mesh
void rewriteBoundary
(
    const bool dryrun,
    const IOobject& io,
    const fileName& regionPrefix,
    HashTable<word>& thisNames,
    HashTable<word>& nbrNames
)
{
    Info<< "Reading boundary from "
        << io.typeFilePath<IOPtrList<entry>>() << nl;

    // Read PtrList of dictionary.
    const word oldTypeName = IOPtrList<entry>::typeName;
    const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
    IOPtrList<entry> patches(io);
    const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
    // Fake type back to what was in field
    const_cast<word&>(patches.type()) = patches.headerClassName();


    // Replace any 'cyclic'
    label nOldCyclics = 0;
    forAll(patches, patchi)
    {
        const dictionary& patchDict = patches[patchi].dict();

        if (patchDict.get<word>("type") == cyclicPolyPatch::typeName)
        {
            if (!patchDict.found("neighbourPatch"))
            {
                Info<< "Patch " << patches[patchi].keyword()
                    << " does not have 'neighbourPatch' entry; assuming it"
                    << " is of the old type." << endl;
                nOldCyclics++;
            }
        }
    }

    Info<< "Detected " << nOldCyclics << " old cyclics." << nl << endl;


    // Save old patches.
    PtrList<entry> oldPatches(patches);

    // Extend
    label nOldPatches = patches.size();
    patches.setSize(nOldPatches+nOldCyclics);

    // Create reordering map
    labelList oldToNew(patches.size());


    // Add new entries
    label addedPatchi = nOldPatches;
    label newPatchi = 0;
    forAll(oldPatches, patchi)
    {
        const dictionary& patchDict = oldPatches[patchi].dict();

        if
        (
            patchDict.get<word>("type") == cyclicPolyPatch::typeName
        )
        {
            const word& name = oldPatches[patchi].keyword();

            if (patchDict.found("neighbourPatch"))
            {
                patches.set(patchi, oldPatches.release(patchi));
                oldToNew[patchi] = newPatchi++;

                // Check if patches come from automatic conversion
                word oldName;

                string::size_type i = name.rfind("_half0");
                if (i != string::npos)
                {
                    oldName = name.substr(0, i);
                    thisNames.insert(oldName, name);
                    Info<< "Detected converted cyclic patch " << name
                        << " ; assuming it originates from " << oldName
                        << endl;
                }
                else
                {
                    i = name.rfind("_half1");
                    if (i != string::npos)
                    {
                        oldName = name.substr(0, i);
                        nbrNames.insert(oldName, name);
                        Info<< "Detected converted cyclic patch " << name
                            << " ; assuming it originates from " << oldName
                            << endl;
                    }
                }
            }
            else
            {
                label nFaces = patchDict.get<label>("nFaces");
                label startFace = patchDict.get<label>("startFace");

                Info<< "Detected old style " << patchDict.get<word>("type")
                    << " patch " << name << " with" << nl
                    << "    nFaces    : " << nFaces << nl
                    << "    startFace : " << startFace << endl;

                word thisName = name + "_half0";
                word nbrName = name + "_half1";

                thisNames.insert(name, thisName);
                nbrNames.insert(name, nbrName);

                // Save current dictionary
                const dictionary patchDict(patches[patchi].dict());

                // Change entry on this side
                patches.set(patchi, oldPatches.release(patchi));
                oldToNew[patchi] = newPatchi++;
                dictionary& thisPatchDict = patches[patchi].dict();
                thisPatchDict.add("neighbourPatch", nbrName);
                thisPatchDict.set("nFaces", nFaces/2);
                patches[patchi].keyword() = thisName;

                // Add entry on other side
                patches.set
                (
                    addedPatchi,
                    new dictionaryEntry
                    (
                        nbrName,
                        dictionary::null,
                        patchDict
                    )
                );
                oldToNew[addedPatchi] = newPatchi++;
                dictionary& nbrPatchDict = patches[addedPatchi].dict();
                nbrPatchDict.set("neighbourPatch", thisName);
                nbrPatchDict.set("nFaces", nFaces/2);
                nbrPatchDict.set("startFace", startFace+nFaces/2);
                patches[addedPatchi].keyword() = nbrName;

                Info<< "Replaced with patches" << nl
                    << patches[patchi].keyword() << " with" << nl
                    << "    nFaces    : "
                    << thisPatchDict.get<label>("nFaces") << nl
                    << "    startFace : "
                    << thisPatchDict.get<label>("startFace") << nl
                    << patches[addedPatchi].keyword() << " with" << nl
                    << "    nFaces    : "
                    << nbrPatchDict.get<label>("nFaces") << nl
                    << "    startFace : "
                    << nbrPatchDict.get<label>("startFace") << nl
                    << endl;

                addedPatchi++;
            }
        }
        else
        {
            patches.set(patchi, oldPatches.release(patchi));
            oldToNew[patchi] = newPatchi++;
        }
    }

    patches.reorder(oldToNew);

    if (returnReduceOr(nOldCyclics))
    {
        if (dryrun)
        {
            //Info<< "-dry-run option: no changes made" << nl << endl;
        }
        else
        {
            if (mvBak(patches.objectPath(), "old"))
            {
                Info<< "Backup to    "
                    << (patches.objectPath() + ".old") << nl;
            }

            Info<< "Write  to    "
                << patches.objectPath() << nl << endl;
            patches.write();
        }
    }
    else
    {
        Info<< "No changes made to boundary file." << nl << endl;
    }
}


void rewriteField
(
    const bool dryrun,
    const Time& runTime,
    const word& fieldName,
    const HashTable<word>& thisNames,
    const HashTable<word>& nbrNames
)
{
    // Read dictionary. (disable class type checking so we can load
    // field)
    Info<< "Loading field " << fieldName << endl;
    const word oldTypeName = IOdictionary::typeName;
    const_cast<word&>(IOdictionary::typeName) = word::null;

    IOdictionary fieldDict
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    const_cast<word&>(IOdictionary::typeName) = oldTypeName;
    // Fake type back to what was in field
    const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();



    dictionary& boundaryField = fieldDict.subDict("boundaryField");

    bool hasChange = false;

    forAllConstIters(thisNames, iter)
    {
        const word& patchName = iter.key();
        const word& newName = iter.val();

        Info<< "Looking for entry for patch " << patchName << endl;

        // Find old patch name either direct or through wildcards
        // Find new patch name direct only

        if
        (
            boundaryField.found(patchName)
        && !boundaryField.found(newName, keyType::LITERAL)
        )
        {
            Info<< "    Changing entry " << patchName << " to " << newName
                << endl;

            dictionary& patchDict = boundaryField.subDict(patchName);

            if (patchDict.found("value"))
            {
                // Remove any value field since wrong size.
                patchDict.remove("value");
            }


            boundaryField.changeKeyword(patchName, newName);
            boundaryField.add
            (
                nbrNames[patchName],
                patchDict
            );
            Info<< "    Adding entry " << nbrNames[patchName] << endl;

            hasChange = true;
        }
    }

    //Info<< "New boundaryField:" << boundaryField << endl;

    if (returnReduceOr(hasChange))
    {
        if (dryrun)
        {
            //Info<< "-test option: no changes made" << endl;
        }
        else
        {
            if (mvBak(fieldDict.objectPath(), "old"))
            {
                Info<< "Backup to    "
                    << (fieldDict.objectPath() + ".old") << nl;
            }

            Info<< "Write  to    "
                << fieldDict.objectPath() << endl;
            fieldDict.regIOobject::write();
        }
    }
    else
    {
        Info<< "No changes made to field " << fieldName << endl;
    }
    Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Tool to upgrade mesh and fields for split cyclics"
    );

    timeSelector::addOptions();

    argList::addDryRunOption
    (
        "Test only do not change any files"
    );
    argList::addBoolOption
    (
        "enableFunctionEntries",
        "Enable expansion of dictionary directives - #include, #codeStream etc"
    );
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"


    // Make sure we do not use the master-only reading since we read
    // fields (different per processor) as dictionaries.
    IOobject::fileModificationChecking = IOobject::timeStamp;


    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool dryrun = args.dryRun();
    if (dryrun)
    {
        Info<< "-dry-run option: no changes made" << nl << endl;
    }
    const bool enableEntries = args.found("enableFunctionEntries");

    fileName regionPrefix;
    {
        // Specified region or default region
        #include "getRegionOption.H"

        regionPrefix = polyMesh::regionName(regionName);
    }

    // Per cyclic patch the new name for this side and the other side
    HashTable<word> thisNames;
    HashTable<word> nbrNames;

    // Rewrite constant boundary file. Return any patches that have been split.
    IOobject io
    (
        "boundary",
        runTime.constant(),
        polyMesh::meshSubDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (io.typeHeaderOk<IOPtrList<entry>>(false))
    {
        rewriteBoundary
        (
            dryrun,
            io,
            regionPrefix,
            thisNames,
            nbrNames
        );
    }



    // Convert any fields

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time: " << runTime.timeName() << endl;

        // See if mesh in time directory
        IOobject io
        (
            "boundary",
            runTime.timeName(),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (io.typeHeaderOk<IOPtrList<entry>>(false))
        {
            rewriteBoundary
            (
                dryrun,
                io,
                regionPrefix,
                thisNames,
                nbrNames
            );
        }


        IOobjectList objects(runTime, runTime.timeName());


        const int oldFlag = entry::disableFunctionEntries;
        if (!enableEntries)
        {
            // By default disable dictionary expansion for fields
            entry::disableFunctionEntries = 1;
        }

        // Rewriting of fields:

        #undef  rewriteFields
        #define rewriteFields(FieldType)                                      \
        for (const word& fieldName : objects.sortedNames<FieldType>())        \
        {                                                                     \
            rewriteField(dryrun, runTime, fieldName, thisNames, nbrNames);    \
        }

        // volFields
        // ~~~~~~~~~

        rewriteFields(volScalarField);
        rewriteFields(volVectorField);
        rewriteFields(volSphericalTensorField);
        rewriteFields(volSymmTensorField);
        rewriteFields(volTensorField);

        // pointFields
        // ~~~~~~~~~~~

        rewriteFields(pointScalarField);
        rewriteFields(pointVectorField);
        rewriteFields(pointSphericalTensorField);
        rewriteFields(pointSymmTensorField);
        rewriteFields(pointTensorField);


        // surfaceFields
        // ~~~~~~~~~~~~~

        rewriteFields(surfaceScalarField);
        rewriteFields(surfaceVectorField);
        rewriteFields(surfaceSphericalTensorField);
        rewriteFields(surfaceSymmTensorField);
        rewriteFields(surfaceTensorField);

        #undef rewriteFields

        entry::disableFunctionEntries = oldFlag;
    }

    return 0;
}


// ************************************************************************* //
