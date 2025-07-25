/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "mapFields.H"
#include "meshToMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mapFields, 0);
    addToRunTimeSelectionTable(functionObject, mapFields, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mapFields::createInterpolation
(
    const dictionary& dict
)
{
    const fvMesh& meshTarget = mesh_;
    const word mapRegionName(dict.get<word>("mapRegion"));

    Info<< name() << ':' << nl
        << "    Reading mesh " << mapRegionName << endl;

    mapRegionPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                mapRegionName,
                meshTarget.time().constant(),
                meshTarget.time(),
                IOobject::MUST_READ
           )
        )
    );

    const fvMesh& mapRegion = mapRegionPtr_();
    const word mapMethodName(dict.get<word>("mapMethod"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Lookup corresponding AMI method
    word patchMapMethodName = meshToMesh::interpolationMethodAMI(mapMethod);

    // Optionally override
    if (dict.readIfPresent("patchMapMethod", patchMapMethodName))
    {
        Info<< "    Patch mapping method: " << patchMapMethodName << endl;
    }

    Info<< "    Creating mesh to mesh interpolation" << endl;

    if (dict.get<bool>("consistent"))
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName
            )
        );
    }
    else
    {
        HashTable<word> patchMap;

        bool createPatchMap = dict.getOrDefault<bool>("createPatchMap", false);

        const entry* ePtr = dict.findEntry("patchMap");

        if (createPatchMap && ePtr)
        {
            FatalIOErrorInFunction(dict)
                << "Requested 'createPatchMap' but 'patchMap' entry provided. "
                << "Please remove one of these entries"
                << exit(FatalIOError);
        }

        if (!createPatchMap && !ePtr)
        {
            FatalIOErrorInFunction(dict)
                << "Either the 'createPatchMap' or 'patchMap' entry must be "
                << "provided when not using the 'consistent' mapping option"
                << exit(FatalIOError);
        }

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const polyBoundaryMesh& pbmT = mapRegion.boundaryMesh();
        DynamicList<label> cuttingIndices;

        if (createPatchMap)
        {
            Info<< "    Creating patchMap and cuttingPatches" << endl;

            if (dict.found("cuttingPatches"))
            {
                WarningInFunction
                    << "Ignoring user supplied cuttingPatches when "
                    << "createPatchMap option is active"
                    << endl;
            }

            forAll(pbmT, patchiT)
            {
                const polyPatch& ppT = pbmT[patchiT];

                if (!polyPatch::constraintType(ppT.type()))
                {
                    const word& patchNameT = ppT.name();
                    const label patchi = pbm.findPatchID(patchNameT);

                    if (patchi == -1)
                    {
                        Info<< "        Adding to cuttingPatches: "
                            << ppT.name() << endl;

                        cuttingIndices.push_back(ppT.index());
                    }
                    else
                    {
                        if (returnReduce(ppT.size(), sumOp<label>()) > 0)
                        {
                            Info<< "        Adding to patchMap: " << ppT.name()
                                << endl;

                            patchMap.set(patchNameT, patchNameT);
                        }
                    }
                }
            }
        }
        else
        {
            // Read patch map
            patchMap = HashTable<word>(ePtr->stream());

            // Read cutting patches using wordRe's
            const wordRes cuttingPatchRes = dict.get<wordRes>("cuttingPatches");
            cuttingIndices.push_back(pbmT.indices(cuttingPatchRes));
        }

        const wordList cuttingPatches(pbmT.names(), cuttingIndices);

        // Checks
        {
            // Find patch names that exist on target mesh that are not included
            // in the patchMap
            wordHashSet unknownPatchNames;
            for (const auto& ppT : pbmT)
            {
                if
                (
                    !polyPatch::constraintType(ppT.type())
                 && !patchMap.found(ppT.name())
                 && returnReduce(ppT.size(), sumOp<label>()) > 0
                )
                {
                    unknownPatchNames.insert(ppT.name());
                }
            }

            for (const label patchiT : cuttingIndices)
            {
                const word& patchName = pbmT[patchiT].name();

                unknownPatchNames.erase(patchName);

                if (patchMap.found(patchName))
                {
                    Info<< "        Removing cutting patch from patchMap: "
                        << patchName << endl;

                    patchMap.erase(patchName);
                }
            }

            if (unknownPatchNames.size())
            {
                FatalErrorInFunction
                    << "Patches not present in source mesh. "
                    << "Add to cuttingPatches? Patches: "
                    << unknownPatchNames
                    << exit(FatalError);
            }
        }

        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName,
                patchMap,
                cuttingPatches
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::mapFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mapRegionPtr_(),
    interpPtr_(),
    fieldNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mapFields::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readEntry("fields", fieldNames_);
        createInterpolation(dict);

        return true;
    }

    return false;
}


bool Foam::functionObjects::mapFields::execute()
{
    Log << type() << " " << name() << " execute:" << nl;

    bool ok = false;

    ok = mapFieldType<scalar>() || ok;
    ok = mapFieldType<vector>() || ok;
    ok = mapFieldType<sphericalTensor>() || ok;
    ok = mapFieldType<symmTensor>() || ok;
    ok = mapFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }
    return true;
}


bool Foam::functionObjects::mapFields::write()
{
    Log << type() << " " << name() << " write:" << nl;

    bool ok = false;

    ok = writeFieldType<scalar>() || ok;
    ok = writeFieldType<vector>() || ok;
    ok = writeFieldType<sphericalTensor>() || ok;
    ok = writeFieldType<symmTensor>() || ok;
    ok = writeFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }

    return true;
}


// ************************************************************************* //
