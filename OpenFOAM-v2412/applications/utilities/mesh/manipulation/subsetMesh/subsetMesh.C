/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    subsetMesh

Group
    grpMeshManipulationUtilities

Description
    Create a mesh subset for a particular region of interest based on a
    cellSet or cellZone.

    See setSet/topoSet utilities on how to define select cells based on
    various shapes.

    Will subset all points, faces and cells needed to make a sub-mesh,
    but not preserve attached boundary types.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubsetter.H"  // Not fvMeshSubset (need two-step subsetting)
#include "argList.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "ReadFields.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Get the exposed patchId or define the exposedPatchName in fvMeshSubset
label getExposedPatchId(const polyMesh& mesh, const word& patchName)
{
    const label patchId = mesh.boundaryMesh().findPatchID(patchName);

    if (patchId == -1)
    {
        fvMeshSubset::exposedPatchName = patchName;
    }

    Info<< "Adding exposed internal faces to "
        << (patchId == -1 ? "new" : "existing")
        << " patch: " << patchName << nl << endl;

    return patchId;
}


labelList nearestPatch(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Count number of faces in exposedPatchIDs
    label nFaces = 0;
    for (const label patchi : patchIDs)
    {
        nFaces += pbm[patchi].size();
    }

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh.nCells());
    List<topoDistanceData<label>> faceData(mesh.nFaces());

    // Start of changes
    labelList patchFaces(nFaces);
    List<topoDistanceData<label>> patchData(nFaces);
    nFaces = 0;
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData<label>(0, patchi);
            ++nFaces;
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );

    // And extract

    labelList nearest(mesh.nFaces());

    bool haveWarned = false;
    forAll(faceData, faceI)
    {
        if (!faceData[faceI].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some faces, e.g. face " << faceI
                    << " at " << mesh.faceCentres()[faceI] << nl
                    << "Using patch " << patchIDs[0] << " as nearest"
                    << endl;
                haveWarned = true;
            }
            nearest[faceI] = patchIDs[0];
        }
        else
        {
            nearest[faceI] = faceData[faceI].data();
        }
    }

    return nearest;
}


//
// Subset DimensionedField/GeometricField
//
template<class FieldType>
PtrList<FieldType> subsetFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objects
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<FieldType>()
    );

    PtrList<FieldType> subFields(fieldObjects.size());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << FieldType::typeName << " (";
        }
        else
        {
            Info<< ' ';
        }
        Info<< io.name();

        FieldType fld
        (
            IOobject
            (
                io.name(),
                baseMesh.time().timeName(),
                baseMesh,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            baseMesh
        );

        subFields.set(nFields, subsetter.interpolate(fld));
        auto& subField = subFields[nFields];
        ++nFields;

        // Subsetting adds 'subset' prefix - rename to match original.
        subField.rename(io.name());
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }

    return subFields;
}


// Subset point fields
template<class FieldType>
PtrList<FieldType> subsetFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objects,
    const pointMesh& pMesh
)
{
    //const fvMesh& baseMesh = subsetter.baseMesh();

    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<FieldType>()
    );

    PtrList<FieldType> subFields(fieldObjects.size());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << FieldType::typeName << " (";
        }
        else
        {
            Info<< ' ';
        }
        Info<< io.name();

        FieldType fld
        (
            IOobject
            (
                io.name(),
                pMesh.thisDb().time().timeName(),
                pMesh.thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            pMesh
        );

        subFields.set(nFields, subsetter.interpolate(fld));
        auto& subField = subFields[nFields];
        ++nFields;

        // Subsetting adds 'subset' prefix - rename to match original.
        subField.rename(io.name());
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }

    return subFields;
}


template<class TopoSet>
void subsetTopoSets
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const labelList& map,
    const fvMesh& subMesh,
    PtrList<TopoSet>& subSets
)
{
    // Read original sets
    PtrList<TopoSet> sets;
    ReadFields<TopoSet>(objects, sets);

    subSets.resize_null(sets.size());

    forAll(sets, seti)
    {
        const TopoSet& set = sets[seti];

        Info<< "Subsetting " << set.type() << " " << set.name() << endl;

        labelHashSet subset(2*min(set.size(), map.size()));

        // Map the data
        forAll(map, i)
        {
            if (set.found(map[i]))
            {
                subset.insert(i);
            }
        }

        subSets.set
        (
            seti,
            new TopoSet
            (
                subMesh,
                set.name(),
                std::move(subset),
                IOobjectOption::AUTO_WRITE
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a mesh subset for a particular region of interest based on a"
        " cellSet or cellZone(s) specified as the first command argument.\n"
        "See setSet/topoSet utilities on how to select cells based on"
        " various shapes."
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addArgument
    (
        "cell-selection",
        "The cellSet name, but with the -zone option this is interpreted"
        " to be a cellZone selection by name(s) or regex.\n"
        "Eg 'mixer' or '( mixer \"moving.*\" )'"
    );

    argList::addOption
    (
        "patch",
        "name",
        "Add exposed internal faces to specified patch"
        " instead of \"oldInternalFaces\""
    );
    argList::addOption
    (
        "patches",
        "wordRes",
        "Add exposed internal faces to closest of specified patches"
        " instead of \"oldInternalFaces\""
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Exclude single or multiple patches from the -patches selection"
    );
    argList::addBoolOption
    (
        "zone",
        "Subset with cellZone(s) instead of cellSet."
        " The command argument may be a list of words or regexs"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "Specify a time for the resulting mesh"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createNamedMesh.H"
    // Make sure pointMesh gets constructed/read as well
    (void)pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

    // arg[1] = word (cellSet) or wordRes (cellZone)
    // const word selectionName = args[1];

    word meshInstance = mesh.pointsInstance();
    word fieldsInstance = runTime.timeName();

    const bool useCellZone = args.found("zone");
    const bool overwrite = args.found("overwrite");
    const bool specifiedInstance = args.readIfPresent
    (
        "resultTime",
        fieldsInstance
    );
    if (specifiedInstance)
    {
        // Set both mesh and field to this time
        meshInstance = fieldsInstance;
    }


    // Default exposed patch id
    labelList exposedPatchIDs(one{}, -1);

    wordRes includePatches, excludePatches;

    if (!args.readListIfPresent<wordRe>("patches", includePatches))
    {
        if (args.found("patch"))
        {
            includePatches.resize(1);
            includePatches.front() = args.get<word>("patch");
        }
    }
    args.readListIfPresent<wordRe>("exclude-patches", excludePatches);

    if (includePatches.size() == 1 && includePatches.front().isLiteral())
    {
        // Select a single patch - no exclude possible
        exposedPatchIDs.front() =
            getExposedPatchId(mesh, includePatches.front());
    }
    else if (!includePatches.empty())
    {
        // Patches selected (sorted order)
        exposedPatchIDs =
            mesh.boundaryMesh().indices(includePatches, excludePatches);

        // Only retain initial, non-processor patches
        const label nNonProcessor
        (
            mesh.boundaryMesh().nNonProcessor()
        );

        forAll(exposedPatchIDs, i)
        {
            if (exposedPatchIDs[i] > nNonProcessor)
            {
                exposedPatchIDs.resize(i);
                break;
            }
        }

        const wordList allPatchNames(mesh.boundaryMesh().names());

        Info<< "Adding exposed internal faces to nearest of patches:" << nl
            << "    include: " << flatOutput(includePatches) << nl
            << "    exclude: " << flatOutput(excludePatches) << nl
            << nl;

        if (exposedPatchIDs.empty())
        {
            FatalErrorInFunction
                << nl << "No patches matched. Patches: "
                << flatOutput(allPatchNames) << nl
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "Adding exposed internal faces to patch \""
            << fvMeshSubset::exposedPatchName
            << "\" (created if necessary)" << nl
            << nl;
    }


    autoPtr<cellSet> cellSetPtr;

    // arg[1] can be a word (for cellSet) or wordRes (for cellZone)

    wordRes zoneNames;
    if (useCellZone)
    {
        wordRes selectionNames(args.getList<wordRe>(1));
        zoneNames.transfer(selectionNames);

        Info<< "Using cellZone " << flatOutput(zoneNames) << nl << endl;

        if (mesh.cellZones().findIndex(zoneNames) == -1)
        {
            FatalErrorInFunction
                << "No cellZones found: " << flatOutput(zoneNames) << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        const word selectionName = args[1];

        Info<< "Using cellSet " << selectionName << nl << endl;

        cellSetPtr.emplace(mesh, selectionName);
    }


    // Two-step mesh subsetting engine
    fvMeshSubsetter subsetter(mesh);

    {
        bitSet selectedCells =
        (
            cellSetPtr
          ? BitSetOps::create(mesh.nCells(), *cellSetPtr)
          : mesh.cellZones().selection(zoneNames)
        );

        if (exposedPatchIDs.size() == 1)
        {
            // Single patch for exposed faces (syncPar)
            subsetter.reset(selectedCells, exposedPatchIDs.front(), true);
        }
        else
        {
            // The nearest patch per face
            labelList nearestExposedPatch(nearestPatch(mesh, exposedPatchIDs));

            labelList exposedFaces
            (
                subsetter.getExposedFaces(selectedCells, true)  // syncPar
            );

            subsetter.setCellSubset
            (
                selectedCells,
                exposedFaces,
                labelUIndList(nearestExposedPatch, exposedFaces)(),
                true  // syncPar
            );
        }

        FixedList<label, 2> cellCount;
        cellCount[0] = subsetter.subMesh().nCells();
        cellCount[1] = mesh.nCells();
        reduce(cellCount, sumOp<label>());

        Info<< "Subset " << cellCount[0] << " of " << cellCount[1]
            << " cells" << nl << nl;
    }


    IOobjectList objects(mesh, runTime.timeName());


    // Read fields and subset
    #undef  createSubsetFields
    #define createSubsetFields(FieldType, Variable)             \
    PtrList<FieldType> Variable                                 \
    (                                                           \
        subsetFields<FieldType>(subsetter, objects)             \
    );


    // Read vol fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    createSubsetFields(volScalarField, vScalarFlds);
    createSubsetFields(volVectorField, vVectorFlds);
    createSubsetFields(volSphericalTensorField, vSphTensorFlds);
    createSubsetFields(volSymmTensorField, vSymmTensorFlds);
    createSubsetFields(volTensorField, vTensorFlds);

    // Read surface fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    createSubsetFields(surfaceScalarField, sScalarFlds);
    createSubsetFields(surfaceVectorField, sVectorFlds);
    createSubsetFields(surfaceSphericalTensorField, sSphTensorFlds);
    createSubsetFields(surfaceSymmTensorField, sSymmTensorFlds);
    createSubsetFields(surfaceTensorField, sTensorFlds);

    // Read dimensioned fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    createSubsetFields(volScalarField::Internal, dScalarFlds);
    createSubsetFields(volVectorField::Internal, dVectorFlds);
    createSubsetFields(volSphericalTensorField::Internal, dSphTensorFlds);
    createSubsetFields(volSymmTensorField::Internal, dSymmTensorFlds);
    createSubsetFields(volTensorField::Internal, dTensorFlds);


    // Read point fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointMesh& pMesh = pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

    #undef  createSubsetFields
    #define createSubsetFields(FieldType, Variable)             \
    PtrList<FieldType> Variable                                 \
    (                                                           \
        subsetFields<FieldType>(subsetter, objects, pMesh)      \
    );

    createSubsetFields(pointScalarField, pScalarFlds);
    createSubsetFields(pointVectorField, pVectorFlds);
    createSubsetFields(pointSphericalTensorField, pSphTensorFlds);
    createSubsetFields(pointSymmTensorField, pSymmTensorFlds);
    createSubsetFields(pointTensorField, pTensorFlds);

    #undef createSubsetFields


    // Read topoSets and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<cellSet> cellSets;
    PtrList<faceSet> faceSets;
    PtrList<pointSet> pointSets;

    {
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
        if (cellSetPtr)
        {
            objects.remove(*cellSetPtr);
        }
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.cellMap(),
            subsetter.subMesh(),
            cellSets
        );
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.faceMap(),
            subsetter.subMesh(),
            faceSets
        );
        subsetTopoSets
        (
            mesh,
            objects,
            subsetter.pointMap(),
            subsetter.subMesh(),
            pointSets
        );
    }


    // Write mesh and fields to new time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (overwrite || specifiedInstance)
    {
        runTime.setTime(instant(fieldsInstance), 0);
        subsetter.subMesh().setInstance(meshInstance);
        topoSet::setInstance(meshInstance, cellSets);
        topoSet::setInstance(meshInstance, faceSets);
        topoSet::setInstance(meshInstance, pointSets);
    }
    else
    {
        ++runTime;
        subsetter.subMesh().setInstance(runTime.timeName());
    }

    Info<< "Writing subsetted mesh and fields to time " << runTime.timeName()
        << endl;
    subsetter.subMesh().write();
    processorMeshes::removeFiles(subsetter.subMesh());

    auto* subPointMeshPtr =
        subsetter.subMesh().thisDb().findObject<pointMesh>
        (
            pointMesh::typeName
        );
    if (subPointMeshPtr)
    {
        pointMesh& subPointMesh = const_cast<pointMesh&>(*subPointMeshPtr);
        subPointMesh.setInstance(subsetter.subMesh().facesInstance());
        subPointMesh.write();
    }


    // Volume fields
    for (const auto& fld : vScalarFlds)     { fld.write(); }
    for (const auto& fld : vVectorFlds)     { fld.write(); }
    for (const auto& fld : vSphTensorFlds)  { fld.write(); }
    for (const auto& fld : vSymmTensorFlds) { fld.write(); }
    for (const auto& fld : vTensorFlds)     { fld.write(); }

    // Surface fields
    for (const auto& fld : sScalarFlds)     { fld.write(); }
    for (const auto& fld : sVectorFlds)     { fld.write(); }
    for (const auto& fld : sSphTensorFlds)  { fld.write(); }
    for (const auto& fld : sSymmTensorFlds) { fld.write(); }
    for (const auto& fld : sTensorFlds)     { fld.write(); }

    // Dimensioned fields
    for (const auto& fld : dScalarFlds)     { fld.write(); }
    for (const auto& fld : dVectorFlds)     { fld.write(); }
    for (const auto& fld : dSphTensorFlds)  { fld.write(); }
    for (const auto& fld : dSymmTensorFlds) { fld.write(); }
    for (const auto& fld : dTensorFlds)     { fld.write(); }

    // Point fields
    for (const auto& fld : pScalarFlds)     { fld.write(); }
    for (const auto& fld : pVectorFlds)     { fld.write(); }
    for (const auto& fld : pSphTensorFlds)  { fld.write(); }
    for (const auto& fld : pSymmTensorFlds) { fld.write(); }
    for (const auto& fld : pTensorFlds)     { fld.write(); }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
