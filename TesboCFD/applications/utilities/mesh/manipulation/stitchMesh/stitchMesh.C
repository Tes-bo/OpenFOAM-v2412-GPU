/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    stitchMesh

Group
    grpMeshManipulationUtilities

Description
    'Stitches' a mesh.

    Takes a mesh and two patches and merges the faces on the two patches
    (if geometrically possible) so the faces become internal.

    Can do
    - 'perfect' match: faces and points on patches align exactly. Order might
    be different though.
    - 'integral' match: where the surfaces on both patches exactly
    match but the individual faces not
    - 'partial' match: where the non-overlapping part of the surface remains
    in the respective patch.

    Note : Is just a front-end to perfectInterface/slidingInterface.

    Comparable to running a meshModifier of the form
    (if masterPatch is called "M" and slavePatch "S"):

    \verbatim
    couple
    {
        type                    slidingInterface;
        masterFaceZoneName      MSMasterZone
        slaveFaceZoneName       MSSlaveZone
        cutPointZoneName        MSCutPointZone
        cutFaceZoneName         MSCutFaceZone
        masterPatchName         M;
        slavePatchName          S;
        typeOfMatch             partial or integral
    }
    \endverbatim


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "slidingInterface.H"
#include "perfectInterface.H"
#include "IOobjectList.H"
#include "ReadFields.H"
#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Checks whether patch present and non-zero
bool checkPatch(const polyBoundaryMesh& bMesh, const word& name)
{
    const label patchi = bMesh.findPatchID(name);

    if (patchi == -1)
    {
        Info<< "No patch " << name << " in mesh" << nl
            << "Known patches: " << bMesh.names() << endl;

        return false;
    }

    if (bMesh[patchi].empty())
    {
        Info<< "Patch " << name << " has zero size" << nl;

        return false;
    }

    return true;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge the faces on specified patches (if geometrically possible)"
        " so that the\n"
        "faces become internal.\n"
        "This utility can be called without arguments (uses stitchMeshDict)"
        " or with\n"
        "two arguments (master/slave patch names)."
    );

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::addOption("dict", "file", "Alternative stitchMeshDict");

    argList::addBoolOption
    (
        "integral",
        "Couple integral master/slave patches (2 argument mode: default)"
    );
    argList::addBoolOption
    (
        "partial",
        "Couple partially overlapping master/slave patches (2 argument mode)"
    );
    argList::addBoolOption
    (
        "perfect",
        "Couple perfectly aligned master/slave patches (2 argument mode)"
    );
    argList::addBoolOption
    (
        "intermediate",
        "Write intermediate stages, not just the final result"
    );
    argList::addOption
    (
        "toleranceDict",
        "file",
        "Dictionary file with tolerances"
    );

    // The arguments are optional (non-mandatory) when using dictionary mode
    argList::noMandatoryArgs();
    argList::addArgument
    (
        "master",
        "The master patch name (non-dictionary mode)"
    );
    argList::addArgument
    (
        "slave",
        "The slave patch name (non-dictionary mode)"
    );

    #include "setRootCase.H"

    // We now handle checking args and general sanity etc.
    const bool useCommandArgs = (args.size() > 1);

    if (useCommandArgs)
    {
        if (args.found("dict"))
        {
            FatalErrorInFunction
                << "Cannot specify both dictionary and command-line arguments"
                << nl
                << endl;
        }

        // If we have arguments - we require all arguments!
        if (!args.check(true, false))
        {
            FatalError.exit();
        }
    }
    else
    {
        // Carp about inapplicable options

        if (args.found("integral"))
        {
            FatalErrorInFunction
                << "Only specify -integral with command-line arguments"
                << endl;
        }

        if (args.found("partial"))
        {
            FatalErrorInFunction
                << "Only specify -partial with command-line arguments"
                << endl;
        }

        if (args.found("perfect"))
        {
            FatalErrorInFunction
                << "Only specify -perfect with command-line arguments"
                << endl;
        }
    }

    #include "createTime.H"
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const bool intermediate = args.found("intermediate");
    const bool overwrite = args.found("overwrite");

    const word dictName("stitchMeshDict");

    // A validated input dictionary
    dictionary validatedDict;

    if (useCommandArgs)
    {
        // Command argument driven:
        const int integralCover = args.found("integral");
        const int partialCover  = args.found("partial");
        const int perfectCover  = args.found("perfect");

        if ((integralCover + partialCover + perfectCover) > 1)
        {
            FatalErrorInFunction
                << "Can only specify one of -integral | -partial | -perfect."
                << nl
                << "Use perfect match option if the patches perfectly align"
                << " (both vertex positions and face centres)" << endl
                << exit(FatalError);
        }

        // Patch names
        const word masterPatchName(args[1]);
        const word slavePatchName(args[2]);

        // Patch names
        Info<< "    " <<  masterPatchName
            << " / " << slavePatchName << nl;

        // Bail out if either patch has problems
        if
        (
            !checkPatch(mesh.boundaryMesh(), masterPatchName)
         || !checkPatch(mesh.boundaryMesh(), slavePatchName)
        )
        {
            FatalErrorInFunction
                << "Cannot continue"
                << exit(FatalError);

            return 1;
        }

        // Input was validated
        dictionary dict;

        if (perfectCover)
        {
            dict.add("match", word("perfect"));
        }
        else if (partialCover)
        {
            dict.add
            (
                "match",
                slidingInterface::typeOfMatchNames[slidingInterface::PARTIAL]
            );
        }
        else
        {
            dict.add
            (
                "match",
                slidingInterface::typeOfMatchNames[slidingInterface::INTEGRAL]
            );
        }

        // Patch names
        dict.add("master", masterPatchName);
        dict.add("slave", slavePatchName);

        validatedDict.add("stitchMesh", dict);
    }
    else
    {
        // dictionary-driven:

        #include "setSystemRunTimeDictionaryIO.H"

        Info<< "Reading " << dictIO.name() << flush;

        IOdictionary stitchDict(dictIO);

        Info<< " with " << stitchDict.size() << " entries" << nl;

        // Suppress duplicate names
        wordHashSet requestedPatches;

        for (const entry& dEntry : stitchDict)
        {
            if (!dEntry.isDict())
            {
                Info<< "Ignoring non-dictionary entry: "
                    << dEntry.keyword() << nl;
                continue;
            }

            const keyType& key = dEntry.keyword();
            const dictionary& dict = dEntry.dict();

            // Match type
            word matchName;
            if (dict.readIfPresent("match", matchName))
            {
                if
                (
                    matchName != "perfect"
                 && !slidingInterface::typeOfMatchNames.found(matchName)
                )
                {
                    Info<< "Error: unknown match type - " << matchName
                        << " should be one of "
                        << slidingInterface::typeOfMatchNames.toc() << nl;
                    continue;
                }
            }

            // Patch names
            const word masterPatchName(dict.get<word>("master"));
            const word slavePatchName(dict.get<word>("slave"));

            // Patch names
            Info<< "    " <<  masterPatchName
                << " / " << slavePatchName << nl;

            if (!requestedPatches.insert(masterPatchName))
            {
                Info<< "Error: patch specified multiple times - "
                    << masterPatchName << nl;
                continue;
            }

            if (!requestedPatches.insert(slavePatchName))
            {
                Info<< "Error: patch specified multiple times - "
                    << slavePatchName << nl;

                requestedPatches.erase(masterPatchName);
                continue;
            }

            // Bail out if either patch has problems
            if
            (
                !checkPatch(mesh.boundaryMesh(), masterPatchName)
             || !checkPatch(mesh.boundaryMesh(), slavePatchName)
            )
            {
                requestedPatches.erase(masterPatchName);
                requestedPatches.erase(slavePatchName);
                continue;
            }

            // Input was validated

            validatedDict.add(key, dict);
        }
    }

    const label nActions = validatedDict.size();

    Info<< nl << nActions << " validated actions" << endl;

    if (!nActions)
    {
        Info<<"\nStopping" << nl << endl;
        return 1;
    }


    // ------------------------------------------
    // This is where the real work begins

    // set up the tolerances for the sliding mesh
    dictionary slidingTolerances;
    if (args.found("toleranceDict"))
    {
        IOdictionary toleranceFile
        (
            IOobject
            (
                args["toleranceDict"],
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        slidingTolerances += toleranceFile;
    }


    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    // Read all current fvFields so they will get mapped
    Info<< "Reading all current volfields" << endl;
    PtrList<volScalarField> volScalarFields;
    ReadFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    ReadFields(mesh, objects, volVectorFields);

    PtrList<volSphericalTensorField> volSphericalTensorFields;
    ReadFields(mesh, objects, volSphericalTensorFields);

    PtrList<volSymmTensorField> volSymmTensorFields;
    ReadFields(mesh, objects, volSymmTensorFields);

    PtrList<volTensorField> volTensorFields;
    ReadFields(mesh, objects, volTensorFields);

    //- Uncomment if you want to interpolate surface fields (usually a bad idea)
    //Info<< "Reading all current surfaceFields" << endl;
    //PtrList<surfaceScalarField> surfaceScalarFields;
    //ReadFields(mesh, objects, surfaceScalarFields);
    //
    //PtrList<surfaceVectorField> surfaceVectorFields;
    //ReadFields(mesh, objects, surfaceVectorFields);
    //
    //PtrList<surfaceTensorField> surfaceTensorFields;
    //ReadFields(mesh, objects, surfaceTensorFields);

    // More precision (for points data)
    IOstream::minPrecision(10);

    polyTopoChanger stitcher(mesh, IOobject::NO_READ);

    // Step through the topology changes
    label actioni = 0;
    for (const entry& dEntry : validatedDict)
    {
        const dictionary& dict = dEntry.dict();

        // Match type
        bool perfect = false;
        slidingInterface::typeOfMatch matchType = slidingInterface::PARTIAL;

        word matchName;
        if (dict.readIfPresent("match", matchName))
        {
            if (matchName == "perfect")
            {
                perfect = true;
            }
            else
            {
                matchType = slidingInterface::typeOfMatchNames[matchName];
            }
        }

        // Patch names
        const word masterPatchName(dict.get<word>("master"));
        const word slavePatchName(dict.get<word>("slave"));

        // Zone names
        const word mergePatchName(masterPatchName + slavePatchName);
        const word cutZoneName(mergePatchName + "CutFaceZone");

        Info<< nl << "========================================" << nl;

        // Information messages
        if (perfect)
        {
            Info<< "Coupling PERFECTLY aligned patches "
                << masterPatchName << " / " << slavePatchName << nl << nl
                << "Resulting (internal) faces in faceZone "
                << cutZoneName << nl << nl
                << "The patch vertices and face centres must align within a"
                << " tolerance relative to the minimum edge length on the patch"
                << nl << endl;
        }
        else if (matchType == slidingInterface::INTEGRAL)
        {
            Info<< "Coupling INTEGRALLY matching of patches "
                << masterPatchName << " / " << slavePatchName << nl << nl
                << "Resulting (internal) faces in faceZone "
                << cutZoneName << nl << nl
                << "The overall area covered by both patches should be"
                << " identical!" << endl
                << "If this is not the case use partial"
                << nl << endl;
        }
        else
        {
            Info<< "Coupling PARTIALLY overlapping patches "
                << masterPatchName << " / " << slavePatchName << nl
                << "Resulting internal faces in faceZone "
                << cutZoneName << nl
                << "Uncovered faces remain in their patch"
                << nl << endl;
        }


        // Master/slave patches
        const polyPatch& masterPatch = mesh.boundaryMesh()[masterPatchName];
        const polyPatch& slavePatch = mesh.boundaryMesh()[slavePatchName];

        mesh.pointZones().clearAddressing();
        mesh.faceZones().clearAddressing();
        mesh.cellZones().clearAddressing();

        // Lists of master and slave faces:
        labelList faceIds;

        // Markup master face ids
        faceIds.resize_nocopy(masterPatch.size());
        Foam::identity(faceIds, masterPatch.start());

        stitcher.clear();
        stitcher.setSize(1);

        if (perfect)
        {
            // Add new (empty) zone for resulting internal faces
            mesh.faceZones()
            (
                cutZoneName,
                true // verbose
            ).resetAddressing(std::move(faceIds), false);


            // Add the perfect interface mesh modifier
            stitcher.set
            (
                0,
                new perfectInterface
                (
                    "couple" + Foam::name(actioni),
                    0,
                    stitcher,
                    cutZoneName,
                    masterPatchName,
                    slavePatchName
                )
            );
        }
        else
        {
            mesh.pointZones()
            (
                mergePatchName + "CutPointZone",
                true // verbose
            ) = labelList();

            mesh.faceZones()
            (
                mergePatchName + "Side0Zone",
                true // verbose
            ).resetAddressing(std::move(faceIds), false);

            // Markup slave face ids
            faceIds.resize_nocopy(slavePatch.size());
            Foam::identity(faceIds, slavePatch.start());

            mesh.faceZones()
            (
                mergePatchName + "Side1Zone",
                true // verbose
            ).resetAddressing(std::move(faceIds), false);

            // Add empty zone for cut faces
            mesh.faceZones()
            (
                cutZoneName,
                true // verbose
            ).resetAddressing(labelList(), false);


            // Add the sliding interface mesh modifier
            stitcher.set
            (
                0,
                new slidingInterface
                (
                    "couple" + Foam::name(actioni),
                    0,
                    stitcher,
                    mergePatchName + "Side0Zone",
                    mergePatchName + "Side1Zone",
                    mergePatchName + "CutPointZone",
                    cutZoneName,
                    masterPatchName,
                    slavePatchName,
                    matchType,      // integral or partial
                    true            // couple/decouple mode
                )
            );

            static_cast<slidingInterface&>(stitcher[0]).setTolerances
            (
                slidingTolerances,
                true
            );
        }

        ++actioni;

        // Advance time for intermediate results or only on final
        if (!overwrite && (intermediate || actioni == nActions))
        {
            ++runTime;
        }

        // Execute all polyMeshModifiers
        autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh(true);

        mesh.movePoints(morphMap->preMotionPoints());

        // Write mesh
        if (overwrite)
        {
            mesh.setInstance(oldInstance);
            stitcher.instance() = oldInstance;
        }

        if (intermediate || actioni == nActions)
        {
            Info<< nl << "Writing polyMesh to time "
                << runTime.timeName() << endl;

            // Bypass runTime write (since only writes at writeTime)
            if
            (
                !runTime.objectRegistry::writeObject
                (
                    IOstreamOption
                    (
                        runTime.writeFormat(),
                        runTime.writeCompression()
                    ),
                    true
                )
            )
            {
                FatalErrorInFunction
                    << "Failed writing polyMesh."
                    << exit(FatalError);
            }

            mesh.faceZones().write();
            mesh.pointZones().write();
            mesh.cellZones().write();

            // Write fields
            runTime.write();
        }
    }

    Info<< "\nEnd\n" <<  endl;

    return 0;
}


// ************************************************************************* //
