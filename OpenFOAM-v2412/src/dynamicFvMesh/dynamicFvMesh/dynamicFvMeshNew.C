/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "staticFvMesh.H"
#include "simplifiedDynamicFvMesh.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    // Note: - do not register the dictionary since dynamicFvMeshes themselves
    // do this.
    // - defaultRegion (region0) gets loaded from constant, other ones
    //   get loaded from constant/<regionname>. Normally we'd use
    //   polyMesh::dbDir() but we haven't got a polyMesh yet ...
    IOobject dictHeader
    (
        "dynamicMeshDict",
        io.time().constant(),
        polyMesh::regionName(io.name()),
        io.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictHeader);

        const word modelType(dict.get<word>("dynamicFvMesh"));

        Info<< "Selecting dynamicFvMesh " << modelType << endl;

        io.time().libs().open
        (
            dict,
            "dynamicFvMeshLibs",
            IOobjectConstructorTablePtr_
        );

        if (!IOobjectConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "dynamicFvMesh table is empty"
                << exit(FatalError);
        }

        auto* doInitCtor = doInitConstructorTable(modelType);
        if (doInitCtor)
        {
            DebugInfo
                << "Constructing dynamicFvMesh with explicit initialisation"
                << endl;

            // Two-step constructor
            // 1. Construct mesh, do not initialise
            autoPtr<dynamicFvMesh> meshPtr(doInitCtor(io, false));

            // 2. Initialise parents and itself
            meshPtr().init(true);

            return meshPtr;
        }

        auto* ctorPtr = IOobjectConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "dynamicFvMesh",
                modelType,
                *IOobjectConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<dynamicFvMesh>(ctorPtr(io));
    }

    DebugInfo
        << "Constructing staticFvMesh with explicit initialisation" << endl;

    // 1. Construct mesh, do not initialise
    autoPtr<dynamicFvMesh> meshPtr(new staticFvMesh(io, false));

    // 2. Initialise parents and itself
    meshPtr().init(true);

    return meshPtr;
}


Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New
(
    const argList& args,
    const Time& runTime
)
{
    const IOobject meshIO
    (
        polyMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    );

    if (args.dryRun() || args.found("dry-run-write"))
    {
        Info
            << "Operating in 'dry-run' mode: case will run for 1 time step.  "
            << "All checks assumed OK on a clean exit" << endl;

        FieldBase::allowConstructFromLargerSize = true;

        // Stop after 1 iteration of the simplified mesh

        if (args.found("dry-run-write"))
        {
            // Using saWriteNow triggers function objects execute(), write()
            runTime.stopAt(Foam::Time::saWriteNow);
        }
        else
        {
            // Using saNoWriteNow triggers function objects execute(),
            // but not write()
            runTime.stopAt(Foam::Time::saNoWriteNow);
        }

        functionObject::outputPrefix = "postProcessing-dry-run";

        return simplifiedMeshes::simplifiedDynamicFvMeshBase::New(meshIO);
    }

    return New(meshIO);
}


// ************************************************************************* //
