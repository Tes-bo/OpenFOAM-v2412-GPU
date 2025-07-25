/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    postProcess

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) or on the command-line for the
    selected set of times on the selected set of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "fileFieldSelection.H"
#include "mapPolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void executeFunctionObjects
(
    const argList& args,
    const Time& runTime,
    fvMesh& mesh,
    const wordList& selectedFields,
    functionObjectList& functions,
    bool lastTime
)
{
    Info<< nl << "Reading fields:" << endl;

    // Read objects in constant directory
    IOobjectList constObjects(mesh, runTime.constant());

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // List of stored objects to clear after executing functionObjects
    DynamicList<regIOobject*> storedObjects
    (
        objects.size() + constObjects.size()
    );

    const auto nameMatcher = [&](const word& name) -> bool
    {
        return selectedFields.contains(name);
    };

    // Read GeometricFields

    #undef  ReadFields
    #define ReadFields(FieldType)                                             \
    readFields<FieldType>(mesh, objects, nameMatcher, storedObjects);

    // Read volFields
    ReadFields(volScalarField);
    ReadFields(volVectorField);
    ReadFields(volSphericalTensorField);
    ReadFields(volSymmTensorField);
    ReadFields(volTensorField);

    // Read internal fields
    ReadFields(volScalarField::Internal);
    ReadFields(volVectorField::Internal);
    ReadFields(volSphericalTensorField::Internal);
    ReadFields(volSymmTensorField::Internal);
    ReadFields(volTensorField::Internal);

    // Read surface fields
    ReadFields(surfaceScalarField);
    ReadFields(surfaceVectorField);
    ReadFields(surfaceSphericalTensorField);
    ReadFields(surfaceSymmTensorField);
    ReadFields(surfaceTensorField);


    // Read point fields.
    const pointMesh& pMesh = pointMesh::New(mesh);
    #undef  ReadPointFields
    #define ReadPointFields(FieldType)                                        \
    readFields<FieldType>(pMesh, objects, nameMatcher, storedObjects);

    ReadPointFields(pointScalarField)
    ReadPointFields(pointVectorField);
    ReadPointFields(pointSphericalTensorField);
    ReadPointFields(pointSymmTensorField);
    ReadPointFields(pointTensorField);


    // Read uniform dimensioned fields

    #undef  ReadUniformFields
    #define ReadUniformFields(FieldType)                                      \
    readUniformFields<FieldType>(constObjects, nameMatcher, storedObjects);

    ReadUniformFields(uniformDimensionedScalarField);
    ReadUniformFields(uniformDimensionedVectorField);
    ReadUniformFields(uniformDimensionedSphericalTensorField);
    ReadUniformFields(uniformDimensionedSymmTensorField);
    ReadUniformFields(uniformDimensionedTensorField);

    Info<< nl << "Executing functionObjects" << endl;

    // Execute the functionObjects in post-processing mode
    functions.execute();

    // Execute the functionObject 'end()' function for the last time
    if (lastTime)
    {
        functions.end();
    }

    while (!storedObjects.empty())
    {
        storedObjects.back()->checkOut();
        storedObjects.pop_back();
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Execute the set of functionObjects specified in the selected"
        " dictionary or on the command-line for the"
        " selected set of times on the selected set of fields"
    );

    timeSelector::addOptions();
    #include "addProfilingOption.H"
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "setRootCase.H"

    if (args.found("list"))
    {
        functionObjectList::list();
        return 0;
    }

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    // Initialize the set of selected fields from the command-line options
    functionObjects::fileFieldSelection fields(mesh);
    if (args.found("fields"))
    {
        fields.resetFieldFilters
        (
            HashSet<wordRe>(args.getList<wordRe>("fields"))
        );
    }
    if (args.found("field"))
    {
        fields.resetFieldFilters(args.get<wordRe>("field"));
    }

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionsDict;

    HashSet<wordRe> fieldFilters(fields.filters());

    // Construct functionObjectList
    autoPtr<functionObjectList> functionsPtr
    (
        functionObjectList::New
        (
            args,
            runTime,
            functionsDict,
            fieldFilters // include any additional command-line fields
        )
    );

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        switch (mesh.readUpdate())
        {
            case polyMesh::POINTS_MOVED:
            {
                functionsPtr->movePoints(mesh);
                break;
            }
            case polyMesh::TOPO_CHANGE:
            case polyMesh::TOPO_PATCH_CHANGE:
            {
                mapPolyMesh mpm(mesh);
                functionsPtr->updateMesh(mpm);
                break;
            }
            case polyMesh::UNCHANGED:
            {
                // No additional work
                break;
            }
            default:
            {
                FatalErrorIn(args.executable())
                    << "Unhandled enumeration"
                    << abort(FatalError);
            }
        }

        fields.resetFieldFilters(fieldFilters);

        fields.updateSelection();

        const bool oldThrowingIOErr = FatalIOError.throwing(true);

        try
        {
            executeFunctionObjects
            (
                args,
                runTime,
                mesh,
                fields.selectionNames(),
                functionsPtr(),
                timei == timeDirs.size()-1
            );

            // Report to output (avoid overwriting values from simulation)
            profiling::print(Info);
        }
        catch (const Foam::IOerror& err)
        {
            Warning << err << endl;
        }

        Info<< endl;

        // Restore previous exception throwing state
        FatalIOError.throwing(oldThrowingIOErr);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
