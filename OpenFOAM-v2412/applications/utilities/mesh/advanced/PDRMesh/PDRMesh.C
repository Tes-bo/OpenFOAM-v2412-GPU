/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    PDRMesh

Group
    grpMeshAdvancedUtilities

Description
    Mesh and field preparation utility for PDR type simulations.

    Reads
    - cellSet giving blockedCells
    - faceSets giving blockedFaces and the patch they should go into

    NOTE: To avoid exposing wrong fields values faceSets should include
    faces contained in the blockedCells cellset.

    - coupledFaces reads coupledFacesSet to introduces mixed-coupled
      duplicate baffles

    Subsets out the blocked cells and splits the blockedFaces and updates
    fields.

    The face splitting is done by duplicating the faces. No duplication of
    points for now (so checkMesh will show a lot of 'duplicate face' messages)

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "cellSet.H"
#include "BitOps.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "mapPolyMesh.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "regionSplit.H"
#include "Tuple2.H"
#include "cyclicFvPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void modifyOrAddFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const bool flipFaceFlux,
    const label newPatchi,
    const label zoneID,
    const bool zoneFlip,

    bitSet& modifiedFace
)
{
    if (modifiedFace.set(facei))
    {
        // First usage of face. Modify.
        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                facei,                      // label of face
                own,                        // owner
                -1,                         // neighbour
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
    else
    {
        // Second or more usage of face. Add.
        meshMod.setAction
        (
            polyAddFace
            (
                f,                          // modified face
                own,                        // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


template<class Type>
PtrList<GeometricField<Type, fvPatchField, volMesh>> subsetVolFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objects,
    const label patchi,
    const Type& exposedValue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const fvMesh& baseMesh = subsetter.baseMesh();

    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<GeoField>()
    );

    PtrList<GeoField> subFields(fieldObjects.size());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << GeoField::typeName << " (";
        }
        else
        {
            Info<< ' ';
        }
        Info<< "    " << io.name() << endl;

        // Read unregistered
        IOobject rio(io, IOobjectOption::NO_REGISTER);
        GeoField origField(rio, baseMesh);

        subFields.set(nFields, subsetter.interpolate(origField));
        auto& subField = subFields[nFields];
        ++nFields;


        // Subsetting adds 'subset' prefix. Rename field to be like original.
        subField.rename(io.name());
        subField.writeOpt(IOobjectOption::AUTO_WRITE);


        // Explicitly set exposed faces (in patchi) to exposedValue.
        if (patchi >= 0)
        {
            fvPatchField<Type>& fld = subField.boundaryFieldRef()[patchi];

            const label newStart = fld.patch().patch().start();
            const label oldPatchi = subsetter.patchMap()[patchi];

            if (oldPatchi == -1)
            {
                // New patch. Reset whole value.
                fld = exposedValue;
            }
            else
            {
                // Reset faces that originate from different patch
                // or internal faces.

                const fvPatchField<Type>& origPfld =
                    origField.boundaryField()[oldPatchi];

                const label oldSize = origPfld.size();
                const label oldStart = origPfld.patch().patch().start();

                forAll(fld, j)
                {
                    const label oldFacei = subsetter.faceMap()[newStart+j];

                    if (oldFacei < oldStart || oldFacei >= oldStart+oldSize)
                    {
                        fld[j] = exposedValue;
                    }
                }
            }
        }
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }

    return subFields;
}


template<class Type>
PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>> subsetSurfaceFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objects,
    const label patchi,
    const Type& exposedValue
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> GeoField;

    const fvMesh& baseMesh = subsetter.baseMesh();

    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<GeoField>()
    );

    PtrList<GeoField> subFields(fieldObjects.size());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "Subsetting " << GeoField::typeName << " (";
        }
        else
        {
            Info<< ' ';
        }
        Info<< io.name();

        // Read unregistered
        IOobject rio(io, IOobjectOption::NO_REGISTER);
        GeoField origField(rio, baseMesh);

        subFields.set(nFields, subsetter.interpolate(origField));
        auto& subField = subFields[nFields];
        ++nFields;

        // Subsetting adds 'subset' prefix. Rename field to be like original.
        subField.rename(io.name());
        subField.writeOpt(IOobjectOption::AUTO_WRITE);


        // Explicitly set exposed faces (in patchi) to exposedValue.
        if (patchi >= 0)
        {
            fvsPatchField<Type>& fld = subField.boundaryFieldRef()[patchi];

            const label newStart = fld.patch().patch().start();
            const label oldPatchi = subsetter.patchMap()[patchi];

            if (oldPatchi == -1)
            {
                // New patch. Reset whole value.
                fld = exposedValue;
            }
            else
            {
                // Reset faces that originate from different patch
                // or internal faces.

                const fvsPatchField<Type>& origPfld =
                    origField.boundaryField()[oldPatchi];

                const label oldSize = origPfld.size();
                const label oldStart = origPfld.patch().patch().start();

                forAll(fld, j)
                {
                    const label oldFacei = subsetter.faceMap()[newStart+j];

                    if (oldFacei < oldStart || oldFacei >= oldStart+oldSize)
                    {
                        fld[j] = exposedValue;
                    }
                }
            }
        }
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }

    return subFields;
}


// Faces introduced into zero-sized patches don't get a value at all.
// This is hack to set them to an initial value.
template<class GeoField>
void initCreatedPatches
(
    const fvMesh& mesh,
    const mapPolyMesh& map,
    const typename GeoField::value_type initValue
)
{
    for (const GeoField& field : mesh.objectRegistry::csorted<GeoField>())
    {
        auto& fieldBf = const_cast<GeoField&>(field).boundaryFieldRef();

        forAll(fieldBf, patchi)
        {
            if (map.oldPatchSizes()[patchi] == 0)
            {
                // Not mapped.
                fieldBf[patchi] = initValue;

                if (fieldBf[patchi].fixesValue())
                {
                    fieldBf[patchi] == initValue;
                }
            }
        }
    }
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

        Info<< "Subsetting " << set.type() << ' ' << set.name() << endl;

        labelHashSet subset(2*min(set.size(), map.size()));

        // Map the data
        forAll(map, i)
        {
            if (set.contains(map[i]))
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


void createCoupledBaffles
(
    fvMesh& mesh,
    const labelList& coupledWantedPatch,
    polyTopoChange& meshMod,
    bitSet& modifiedFace
)
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(coupledWantedPatch, facei)
    {
        if (coupledWantedPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];
            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            // Use owner side of face
            modifyOrAddFace
            (
                meshMod,
                f,                          // modified face
                facei,                      // label of face
                mesh.faceOwner()[facei],    // owner
                false,                      // face flip
                coupledWantedPatch[facei],  // patch for face
                zoneID,                     // zone for face
                zoneFlip,                   // face flip in zone
                modifiedFace                // modify or add status
            );

            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }
                // Use neighbour side of face
                modifyOrAddFace
                (
                    meshMod,
                    f.reverseFace(),            // modified face
                    facei,                      // label of face
                    mesh.faceNeighbour()[facei],// owner
                    false,                      // face flip
                    coupledWantedPatch[facei],  // patch for face
                    zoneID,                     // zone for face
                    zoneFlip,                   // face flip in zone
                    modifiedFace                // modify or add status
                );
            }
        }
    }
}


void createCyclicCoupledBaffles
(
    fvMesh& mesh,
    const labelList& cyclicMasterPatch,
    const labelList& cyclicSlavePatch,
    polyTopoChange& meshMod,
    bitSet& modifiedFace
)
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(cyclicMasterPatch, facei)
    {
        if (cyclicMasterPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];

            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            modifyOrAddFace
            (
                meshMod,
                f.reverseFace(),                    // modified face
                facei,                              // label of face
                mesh.faceNeighbour()[facei],        // owner
                false,                              // face flip
                cyclicMasterPatch[facei],           // patch for face
                zoneID,                             // zone for face
                zoneFlip,                           // face flip in zone
                modifiedFace                        // modify or add
            );
        }
    }

    forAll(cyclicSlavePatch, facei)
    {
        if (cyclicSlavePatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];
            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }
            // Use owner side of face
                modifyOrAddFace
                (
                    meshMod,
                    f,                          // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],    // owner
                    false,                      // face flip
                    cyclicSlavePatch[facei],    // patch for face
                    zoneID,                     // zone for face
                    zoneFlip,                   // face flip in zone
                    modifiedFace                // modify or add status
                );
            }
        }
    }
}


void createBaffles
(
    fvMesh& mesh,
    const labelList& wantedPatch,
    polyTopoChange& meshMod
)
{
    const faceZoneMesh& faceZones = mesh.faceZones();
    Info << "faceZone:createBaffle " << faceZones << endl;
    forAll(wantedPatch, facei)
    {
        if (wantedPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];

            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                          // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],    // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    wantedPatch[facei],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );

            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }

                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        mesh.faceNeighbour()[facei],// owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        facei,                      // masterFaceID,
                        false,                      // face flip
                        wantedPatch[facei],         // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }
    }
}


// Wrapper around find patch. Also makes sure same patch in parallel.
label findPatch(const polyBoundaryMesh& patches, const word& patchName)
{
    const label patchi = patches.findPatchID(patchName);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Illegal patch " << patchName
            << nl << "Valid patches are " << patches.names()
            << exit(FatalError);
    }

    // Check same patch for all procs
    {
        const label newPatchi = returnReduce(patchi, minOp<label>());

        if (newPatchi != patchi)
        {
            FatalErrorInFunction
                << "Patch " << patchName
                << " should have the same patch index on all processors." << nl
                << "On my processor it has index " << patchi
                << " ; on some other processor it has index " << newPatchi
                << exit(FatalError);
        }
    }
    return patchi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mesh and field preparation utility for PDR type simulations."
    );
    #include "addOverwriteOption.H"

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~

    IOdictionary dict
    (
        IOobject
        (
            "PDRMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Per faceSet the patch to put the baffles into
    const List<Pair<word>> setsAndPatches(dict.lookup("blockedFaces"));

    // Per faceSet the patch to put the coupled baffles into
    DynamicList<FixedList<word, 3>> coupledAndPatches(10);

    const dictionary& functionDicts = dict.subDict("coupledFaces");

    for (const entry& dEntry : functionDicts)
    {
        if (!dEntry.isDict())  // Safety
        {
            continue;
        }

        const word& key = dEntry.keyword();
        const dictionary& dict = dEntry.dict();

        const word cyclicName = dict.get<word>("cyclicMasterPatch");
        const word wallName = dict.get<word>("wallPatch");
        FixedList<word, 3> nameAndType;
        nameAndType[0] = key;
        nameAndType[1] = wallName;
        nameAndType[2] = cyclicName;
        coupledAndPatches.append(nameAndType);
    }

    forAll(setsAndPatches, setI)
    {
        Info<< "Faces in faceSet " << setsAndPatches[setI][0]
            << " become baffles in patch " << setsAndPatches[setI][1]
            << endl;
    }

    forAll(coupledAndPatches, setI)
    {
        Info<< "Faces in faceSet " << coupledAndPatches[setI][0]
            << " become coupled baffles in patch " << coupledAndPatches[setI][1]
            << endl;
    }

    // All exposed faces that are not explicitly marked to be put into a patch
    const word defaultPatch(dict.get<word>("defaultPatch"));

    Info<< "Faces that get exposed become boundary faces in patch "
        << defaultPatch << endl;

    const word blockedSetName(dict.get<word>("blockedCells"));

    Info<< "Reading blocked cells from cellSet " << blockedSetName
        << endl;

    const bool overwrite = args.found("overwrite");


    // Read faceSets, lookup patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Check that face sets don't have coincident faces
    labelList wantedPatch(mesh.nFaces(), -1);
    forAll(setsAndPatches, setI)
    {
        faceSet fSet(mesh, setsAndPatches[setI][0]);

        label patchi = findPatch
        (
            mesh.boundaryMesh(),
            setsAndPatches[setI][1]
        );

        for (const label facei : fSet)
        {
            if (wantedPatch[facei] != -1)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " is in faceSet " << setsAndPatches[setI][0]
                    << " destined for patch " << setsAndPatches[setI][1]
                    << " but also in patch " << wantedPatch[facei]
                    << exit(FatalError);
            }
            wantedPatch[facei] = patchi;
        }
    }

    // Per face the patch for coupled baffle or -1.
    labelList coupledWantedPatch(mesh.nFaces(), -1);
    labelList cyclicWantedPatch_half0(mesh.nFaces(), -1);
    labelList cyclicWantedPatch_half1(mesh.nFaces(), -1);

    forAll(coupledAndPatches, setI)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const label cyclicId =
            findPatch(patches, coupledAndPatches[setI][2]);

        const label cyclicSlaveId = findPatch
        (
            patches,
            refCast<const cyclicFvPatch>
            (
                mesh.boundary()[cyclicId]
            ).neighbFvPatch().name()
        );

        faceSet fSet(mesh, coupledAndPatches[setI][0]);
        label patchi = findPatch(patches, coupledAndPatches[setI][1]);

        for (const label facei : fSet)
        {
            if (coupledWantedPatch[facei] != -1)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " is in faceSet " << coupledAndPatches[setI][0]
                    << " destined for patch " << coupledAndPatches[setI][1]
                    << " but also in patch " << coupledWantedPatch[facei]
                    << exit(FatalError);
            }

            coupledWantedPatch[facei] = patchi;
            cyclicWantedPatch_half0[facei] = cyclicId;
            cyclicWantedPatch_half1[facei] = cyclicSlaveId;
        }
    }

    // Exposed faces patch
    label defaultPatchi = findPatch(mesh.boundaryMesh(), defaultPatch);


    //
    // Removing blockedCells (like subsetMesh)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //

    // Mesh subsetting engine
    fvMeshSubset subsetter
    (
        mesh,
        BitSetOps::create
        (
            mesh.nCells(),
            cellSet(mesh, blockedSetName), // Blocked cells as labelHashSet
            false  // on=false: invert logic => retain the unblocked cells
        ),
        defaultPatchi,
        true
    );


    // Subset wantedPatch. Note that might also include boundary faces
    // that have been exposed by subsetting.
    wantedPatch = IndirectList<label>(wantedPatch, subsetter.faceMap())();

    coupledWantedPatch = IndirectList<label>
    (
        coupledWantedPatch,
        subsetter.faceMap()
    )();

    cyclicWantedPatch_half0 = IndirectList<label>
    (
        cyclicWantedPatch_half0,
        subsetter.faceMap()
    )();

    cyclicWantedPatch_half1 = IndirectList<label>
    (
        cyclicWantedPatch_half1,
        subsetter.faceMap()
    )();

    // Read all fields in time and constant directories
    IOobjectList objects(mesh, runTime.timeName());
    {
        IOobjectList timeObjects(mesh, mesh.facesInstance());

        // Transfer specific types
        forAllIters(timeObjects, iter)
        {
            autoPtr<IOobject> objPtr(timeObjects.remove(iter));
            const auto& obj = *objPtr;

            if
            (
                obj.isHeaderClass<volScalarField>()
             || obj.isHeaderClass<volVectorField>()
             || obj.isHeaderClass<volSphericalTensorField>()
             || obj.isHeaderClass<volTensorField>()
             || obj.isHeaderClass<volSymmTensorField>()
             || obj.isHeaderClass<surfaceScalarField>()
             || obj.isHeaderClass<surfaceVectorField>()
             || obj.isHeaderClass<surfaceSphericalTensorField>()
             || obj.isHeaderClass<surfaceSymmTensorField>()
             || obj.isHeaderClass<surfaceTensorField>()
            )
            {
                objects.add(objPtr);
            }
        }
    }

    // Read vol fields and subset.
    PtrList<volScalarField> scalarFlds
    (
        subsetVolFields<scalar>
        (
            subsetter,
            objects,
            defaultPatchi,
            scalar(Zero)
        )
    );

    PtrList<volVectorField> vectorFlds
    (
        subsetVolFields<vector>
        (
            subsetter,
            objects,
            defaultPatchi,
            vector(Zero)
        )
    );

    PtrList<volSphericalTensorField> sphTensorFlds
    (
        subsetVolFields<sphericalTensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            sphericalTensor(Zero)
        )
    );

    PtrList<volSymmTensorField> symmTensorFlds
    (
        subsetVolFields<symmTensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            symmTensor(Zero)
        )
    );

    PtrList<volTensorField> tensorFlds
    (
        subsetVolFields<tensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            tensor(Zero)
        )
    );

    // Read surface fields and subset.
    PtrList<surfaceScalarField> surfScalarFlds
    (
        subsetSurfaceFields<scalar>
        (
            subsetter,
            objects,
            defaultPatchi,
            scalar(Zero)
        )
    );

    PtrList<surfaceVectorField> surfVectorFlds
    (
        subsetSurfaceFields<vector>
        (
            subsetter,
            objects,
            defaultPatchi,
            vector(Zero)
        )
    );

    PtrList<surfaceSphericalTensorField> surfSphericalTensorFlds
    (
        subsetSurfaceFields<sphericalTensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            sphericalTensor(Zero)
        )
    );

    PtrList<surfaceSymmTensorField> surfSymmTensorFlds
    (
        subsetSurfaceFields<symmTensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            symmTensor(Zero)
        )
    );

    PtrList<surfaceTensorField> surfTensorFlds
    (
        subsetSurfaceFields<tensor>
        (
            subsetter,
            objects,
            defaultPatchi,
            tensor(Zero)
        )
    );


    // Set handling
    PtrList<cellSet> cellSets;
    PtrList<faceSet> faceSets;
    PtrList<pointSet> pointSets;
    {
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
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


    if (!overwrite)
    {
        ++runTime;
    }

    Info<< "Writing mesh without blockedCells to time "
        << runTime.value() << endl;

    subsetter.subMesh().write();


    //
    // Splitting blockedFaces
    // ~~~~~~~~~~~~~~~~~~~~~~
    //

    // Synchronize wantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        wantedPatch,
        maxEqOp<label>()
    );

    // Synchronize coupledWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        coupledWantedPatch,
        maxEqOp<label>()
    );

    // Synchronize cyclicWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half0,
        maxEqOp<label>()
    );

    // Synchronize cyclicWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half1,
        maxEqOp<label>()
    );

    // Topochange container
    polyTopoChange meshMod(subsetter.subMesh());


    // Whether first use of face (modify) or consecutive (add)
    bitSet modifiedFace(mesh.nFaces());

    // Create coupled wall-side baffles
    createCoupledBaffles
    (
        subsetter.subMesh(),
        coupledWantedPatch,
        meshMod,
        modifiedFace
    );

    // Create coupled master/slave cyclic baffles
    createCyclicCoupledBaffles
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half0,
        cyclicWantedPatch_half1,
        meshMod,
        modifiedFace
    );

    // Create wall baffles
    createBaffles
    (
        subsetter.subMesh(),
        wantedPatch,
        meshMod
    );

    if (!overwrite)
    {
        ++runTime;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> mapPtr =
        meshMod.changeMesh(subsetter.subMesh(), false);
    mapPolyMesh& map = *mapPtr;

    // Update fields
    subsetter.subMesh().updateMesh(map);

    // Fix faces that get mapped to zero-sized patches (these don't get any
    // value)
    initCreatedPatches<volScalarField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volVectorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volSphericalTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volSymmTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );

    initCreatedPatches<surfaceScalarField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceVectorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceSphericalTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceSymmTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );

    // Update numbering of topoSets
    topoSet::updateMesh(subsetter.subMesh().facesInstance(), map, cellSets);
    topoSet::updateMesh(subsetter.subMesh().facesInstance(), map, faceSets);
    topoSet::updateMesh(subsetter.subMesh().facesInstance(), map, pointSets);


    // Move mesh (since morphing might not do this)
    if (map.hasMotionPoints())
    {
        subsetter.subMesh().movePoints(map.preMotionPoints());
    }

    Info<< "Writing mesh with split blockedFaces to time " << runTime.value()
        << endl;

    subsetter.subMesh().write();


    //
    // Removing inaccessible regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //

    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.
    regionSplit cellRegion(subsetter.subMesh());

    if (cellRegion.nRegions() > 1)
    {
        WarningInFunction
            << "Removing blocked faces and cells created "
            << cellRegion.nRegions()
            << " regions that are not connected via a face." << nl
            << "    This is not supported in solvers." << nl
            << "    Use" << nl << nl
            << "    splitMeshRegions <root> <case> -largestOnly" << nl << nl
            << "    to extract a single region of the mesh." << nl
            << "    This mesh will be written to a new timedirectory"
            << " so might have to be moved back to constant/" << nl
            << endl;

        const word startFrom(runTime.controlDict().get<word>("startFrom"));

        if (startFrom != "latestTime")
        {
            WarningInFunction
                << "To run splitMeshRegions please set your"
                << " startFrom entry to latestTime" << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
