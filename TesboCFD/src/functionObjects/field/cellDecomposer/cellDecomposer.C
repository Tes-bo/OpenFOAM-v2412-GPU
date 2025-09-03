/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "cellDecomposer.H"
#include "addToRunTimeSelectionTable.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "tetDecomposer.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "ReadFields.H"
#include "surfaceFields.H"
#include "PackedBoolList.H"
#include "fvMeshTools.H"
#include "cellSetOption.H"
#include "cellBitSet.H"
#include "cellSet.H"
#include "hexMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetMatcher.H"

#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cellDecomposer, 0);
    addToRunTimeSelectionTable(functionObject, cellDecomposer, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cellDecomposer::makeMesh
(
    const dictionary& dict,
    const word& regionName
)
{
    Info<< name() << ':' << nl
        << "    Decomposing cells to region " << regionName << endl;

    Info<< incrIndent;

    const fv::cellSetOption::selectionModeType selectionMode
    (
        fv::cellSetOption::selectionModeTypeNames_.get
        (
            "selectionMode",
            dict
        )
    );


    // Start off from existing mesh
    polyTopoChange meshMod(mesh_);

    tetDecompPtr_.reset(new Foam::tetDecomposer(mesh_));

    // Default values
    bitSet decomposeCell(mesh_.nCells());
    autoPtr<bitSet> decomposeFacePtr;
    tetDecomposer::decompositionType decompType
    (
        tetDecomposer::FACE_CENTRE_TRIS
    );


    switch (selectionMode)
    {
        case fv::cellSetOption::smAll:
        {
            Info<< indent << "- selecting all cells" << endl;

            decomposeCell = true;
            break;
        }
        case fv::cellSetOption::smGeometric:
        {
            Info<< indent << "- selecting cells geometrically" << endl;

            decomposeCell =
                cellBitSet::select(mesh_, dict.subDict("selection"), true);
            break;
        }
        case fv::cellSetOption::smCellSet:
        {
            const word selectionName(dict.get<word>("cellSet"));

            Info<< indent
                << "- selecting cells using cellSet "
                << selectionName << endl;

            decomposeCell.set
            (
                cellSet(mesh_, selectionName, IOobject::MUST_READ).toc()
            );
            break;
        }
        case fv::cellSetOption::smCellZone:
        {
            wordRes selectionNames;
            if
            (
                !dict.readIfPresent("cellZones", selectionNames)
             || selectionNames.empty()
            )
            {
                selectionNames.resize(1);
                dict.readEntry("cellZone", selectionNames.first());
            }

            Info<< indent
                << "- selecting cells using cellZones "
                << flatOutput(selectionNames) << nl;

            const auto& zones = mesh_.cellZones();

            // Also handles groups, multiple zones etc ...
            labelList zoneIDs = zones.indices(selectionNames);

            if (zoneIDs.empty())
            {
                FatalErrorInFunction
                    << "No matching cellZones: "
                    << flatOutput(selectionNames) << nl
                    << "Valid zones : "
                    << flatOutput(zones.names()) << nl
                    << "Valid groups: "
                    << flatOutput(zones.groupNames())
                    << nl
                    << exit(FatalError);
            }

            if (zoneIDs.size() == 1)
            {
                decomposeCell.set(zones[zoneIDs.first()]);
                // TBD: Foam::sort(cells_);
            }
            else
            {
                decomposeCell.set(zones.selection(zoneIDs).sortedToc());
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unsupported selectionMode "
                << fv::cellSetOption::selectionModeTypeNames_[selectionMode]
                << ". Valid selectionMode types are "
                << fv::cellSetOption::selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    word decompTypeName;
    if (dict.readIfPresent("decomposeType", decompTypeName))
    {
        if (decompTypeName == "polyhedral")
        {
            // Automatic selection to generate hex/prism/tet only
            //  - subset decomposeCell to exclude hex/prism/tet
            //  - set faces to FACE_DIAG_QUADS
            // Usually start with cellSetOption::smAll
            decomposeFacePtr.reset(new bitSet(mesh_.nFaces()));
            auto& decomposeFace = decomposeFacePtr();

            decompType = tetDecomposer::FACE_DIAG_QUADS;
            bitSet oldDecomposeCell(decomposeCell);
            decomposeCell = false;

            // Construct shape recognizers
            prismMatcher prism;

            for (const label celli : oldDecomposeCell)
            {
                if
                (
                   !hexMatcher::test(mesh_, celli)
                && !tetMatcher::test(mesh_, celli)
                && !pyrMatcher::test(mesh_, celli)
                && !prism.isA(mesh_, celli)
                )
                {
                    decomposeCell.set(celli);
                    decomposeFace.set(mesh_.cells()[celli]);
                }
            }

            // Sync boundary info
            syncTools::syncFaceList
            (
                mesh_,
                decomposeFace,
                orEqOp<unsigned int>()
            );
        }
        else
        {
            decompType = tetDecomposer::decompositionTypeNames[decompTypeName];
        }
    }


    // Insert all changes to create tets
    if (decomposeFacePtr)
    {
        if (debug)
        {
            OBJstream os(mesh_.time().path()/"orig_faces.obj");
            os.write
            (
                UIndirectList<face>
                (
                    mesh_.faces(),
                    decomposeFacePtr().sortedToc()
                )(),
                mesh_.points(),
                false
            );
            Pout<< "Written " << meshMod.faces().size()
                << " faces to " << os.name() << endl;
        }

        tetDecompPtr_().setRefinement
        (
            decompType, //Foam::tetDecomposer::FACE_CENTRE_TRIS,
            decomposeCell,
            decomposeFacePtr(),
            meshMod
        );
    }
    else
    {
        tetDecompPtr_().setRefinement
        (
            decompType, //Foam::tetDecomposer::FACE_CENTRE_TRIS,
            decomposeCell,
            meshMod
        );
    }


    if (debug)
    {
        OBJstream os(mesh_.time().path()/"faces.obj");
        os.write(meshMod.faces(), pointField(meshMod.points()), false);
        Pout<< "Written " << meshMod.faces().size()
            << " faces to " << os.name() << endl;
    }


    autoPtr<fvMesh> tetMeshPtr;

    mapPtr_ = meshMod.makeMesh
    (
        tetMeshPtr,
        IOobject
        (
            regionName,
            mesh_.facesInstance(),      // same instance as original mesh
            mesh_.time(),               //? why same database? TBD.
            IOobject::READ_IF_PRESENT,  // Read fv* if present
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    );


    //- Update numbering on tet-decomposition engine
    tetDecompPtr_().updateMesh(mapPtr_());

    Info<< indent << "Writing decomposed mesh to "
        << tetMeshPtr().objectRegistry::objectRelPath()
        << endl;
    tetMeshPtr().write();

    fvMeshTools::createDummyFvMeshFiles(mesh_.time(), regionName, true);

    // Store new mesh on object registry
    tetMeshPtr.ptr()->polyMesh::store();

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cellDecomposer::cellDecomposer
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cellDecomposer::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        // Generate new tet equivalent mesh. TBD: why meshObject?
        //meshObjects::tetDecomposition::New(mesh, dict);

        dict_ = dict.optionalSubDict(typeName + "Coeffs");
        dict_.readEntry("mapRegion", mapRegion_);
        dict_.readEntry("fields", fieldNames_);
        makeMesh(dict_, mapRegion_);
    }

    return true;
}


bool Foam::functionObjects::cellDecomposer::execute()
{
    Log << type() << " " << name() << " execute:" << nl;

    bool ok = false;

    if (mesh_.changing())
    {
        // Re-do mesh. Currently rebuilds whole mesh. Could move points only
        // for mesh motion.
        tetDecompPtr_.clear();
        mapPtr_.clear();
        const_cast<Time&>(this->mesh_.time()).erase(mapRegion_);
        makeMesh(dict_, mapRegion_);
    }


    // Look up
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


bool Foam::functionObjects::cellDecomposer::write()
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
