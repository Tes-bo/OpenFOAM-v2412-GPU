/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "columnFvMesh.H"
#include "polyBoundaryMeshEntries.H"
#include "IOobjectList.H"
#include "fieldDictionary.H"
#include "vectorIOField.H"
#include "emptyPolyPatch.H"
#include "topoSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedMeshes
{
    defineTypeNameAndDebug(columnFvMeshInfo, 0);
    defineTypeNameAndDebug(columnFvMesh, 0);

    addToRunTimeSelectionTable
    (
        simplifiedFvMesh,
        columnFvMesh,
        time
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::simplifiedMeshes::columnFvMeshInfo::setPatchEntries
(
    const Time& runTime
)
{
    IOobject boundaryIO
    (
        "boundary",
        localInstance_,
        polyMesh::meshDir(regionName_),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    const wordHashSet constraintPatches(polyPatch::constraintTypes());

    if (boundaryIO.typeHeaderOk<polyBoundaryMesh>(true))
    {
        polyBoundaryMeshEntries allPatchEntries(boundaryIO);

        Info<< "Creating simplified mesh using " << allPatchEntries.path()
            << endl;

        for (const entry& e : allPatchEntries)
        {
            const word type(e.dict().get<word>("type"));

            if (!constraintPatches.found(type))
            {
                if (e.dict().get<label>("nFaces"))
                {
                    ++nPatchWithFace_;
                }
                patchEntries_.add(e.keyword(), e.dict());
            }
        }

        return true;
    }
    else
    {
        // No boundary file - try reading from a field
        IOobjectList objects
        (
            runTime,
            runTime.timeName(),
            polyMesh::regionName(regionName_)
        );

        if (objects.empty())
        {
            // No fields - cannot create a mesh
            return false;
        }

        const IOobject& io = *objects.begin()();

        if (io.instance() == runTime.constant())
        {
            FatalErrorInFunction
                << "No time directories found for field reading"
                << exit(FatalError);
        }

        const fieldDictionary fieldDict(io, io.headerClassName());

        Info<< "Creating simplified mesh from field "
            << fieldDict.objectPath()
            << endl;

        WarningInFunction
            << "All boundaries will be approximated using wall-type patches. "
            << "This may cause your" << nl
            << "    final case to run differently. "
            << "Create your mesh first for improved performance"
            << endl;

        const dictionary& boundaryFieldDict =
            fieldDict.subDict("boundaryField");

        for (const entry& e : boundaryFieldDict)
        {
            const word type(e.dict().get<word>("type"));

            if (simplifiedFvMesh::fvPatchFieldExists(type))
            {
                if (!constraintPatches.found(type))
                {
                    ++nPatchWithFace_;
                    dictionary simplifiedEntries;
                    simplifiedEntries.add("startFace", 0);
                    simplifiedEntries.add("nFaces", 1);
                    simplifiedEntries.add("type", "wall"); // default to wall

                    patchEntries_.add(e.keyword(), simplifiedEntries);
                }
            }
            else
            {
                Info<< "Ignoring unknown patch type " << type << endl;
            }
        }

        return false;
    }
}


void Foam::simplifiedMeshes::columnFvMeshInfo::initialise(const Time& runTime)
{
    DebugInfo << "Constructing 1-D mesh" << nl << endl;

    // Read patches - filter out proc patches for parallel runs
    createFromMesh_ = setPatchEntries(runTime);

    const label nPatch = patchEntries_.size();

    DebugPout << "Read " << nPatch << " patches" << endl;

    /*
               3 ---- 7
          f5   |\     |\   f2
          |    | 2 ---- 6   \
          |    0 |--- 4 |    \
          |     \|     \|    f3
          f4     1 ---- 5

                f0 ----- f1
    */

    // Create mesh with nPatchWithFace hex cells
    // - hex face 1: set internal face to link to next hex cell
    // - hex face 2,3,4,5: set to boundary condition
    //
    // TODO: Coupled patches:
    // - separated/collocated cyclic - faces 2 and 3
    // - rotational cyclic - convert to separated and warn not handled?
    // - proc patches - faces need to be collocated?

    points1D_.setSize(nPatchWithFace_*4 + 4);
    faces1D_.setSize(nPatchWithFace_*5 + 1);

    owner1D_.setSize(faces1D_.size(), label(-1));
    neighbour1D_.setSize(owner1D_.size(), label(-1));

    vector dx(Zero);
    {
        // Determine the mesh bounds
        IOobject pointsIO
        (
            "points",
            localInstance_,
            polyMesh::meshDir(regionName_),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        scalar dxi = 0;
        scalar dyi = 0;
        scalar dzi = 0;
        point origin(Zero);

        if (pointsIO.typeHeaderOk<vectorIOField>(true))
        {
            const pointIOField meshPoints(pointsIO);

            boundBox meshBb(meshPoints, true);

            Info<< "Mesh bounds: " << meshBb << endl;

            origin = meshBb.min();
            vector span = meshBb.span();
            dxi = span.x()/scalar(nPatchWithFace_);
            dyi = span.y();
            dzi = span.z();
        }
        else
        {
            scalar Lref = GREAT;
            origin = point(-Lref, -Lref, -Lref);
            dxi = 2.0*Lref/scalar(nPatchWithFace_);
            dyi = Lref;
            dzi = Lref;
        }

        dx = vector(dxi, 0, 0);
        const vector dy(0, dyi, 0);
        const vector dz(0, 0, dzi);

        // End face
        points1D_[0] = origin;
        points1D_[1] = origin + dy;
        points1D_[2] = origin + dy + dz;
        points1D_[3] = origin + dz;
    }

    label n = 4;
    for (label i = 1; i <= nPatchWithFace_; ++i)
    {
        const vector idx(i*dx);
        points1D_[i*n] = points1D_[0] + idx;
        points1D_[i*n + 1] = points1D_[1] + idx;
        points1D_[i*n + 2] = points1D_[2] + idx;
        points1D_[i*n + 3] = points1D_[3] + idx;
    }

    if (debug) Pout<< "points:" << points1D_ << endl;

    label facei = 0;

    // Internal faces first
    for (label i = 0; i < nPatchWithFace_ - 1; ++i)
    {
        label o = i*4;
        faces1D_[facei] = face({4 + o, 5 + o, 6 + o, 7 + o});
        owner1D_[facei] = i;
        neighbour1D_[facei] = i + 1;
        ++facei;
    }

    // Boundary conditions
    for (label i = 0; i < nPatchWithFace_; ++i)
    {
        label o = i*4;
        faces1D_[facei] = face({0 + o, 4 + o, 7 + o, 3 + o});
        owner1D_[facei] = i;
        ++facei;

        faces1D_[facei] = face({0 + o, 1 + o, 5 + o, 4 + o});
        owner1D_[facei] = i;
        ++facei;

        faces1D_[facei] = face({1 + o, 2 + o, 6 + o, 5 + o});
        owner1D_[facei] = i;
        ++facei;

        faces1D_[facei] = face({3 + o, 7 + o, 6 + o, 2 + o});
        owner1D_[facei] = i;
        ++facei;
    }
    {
        // End caps
        faces1D_[facei] = face({0, 3, 2, 1});
        owner1D_[facei] = 0;
        ++facei;

        label o = 4*nPatchWithFace_;
        faces1D_[facei] = face({0 + o, 1 + o, 2 + o, 3 + o});
        owner1D_[facei] = nPatchWithFace_ - 1;
        ++facei;
    }

    DebugPout
        << "faces:" << faces1D_ << nl
        << "owner:" << owner1D_ << nl
        << "neighbour:" << neighbour1D_
        << endl;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::simplifiedMeshes::columnFvMeshInfo::addLocalPatches
(
    fvMesh& mesh
) const
{
    const label nPatch = patchEntries_.size();

    List<polyPatch*> patches(nPatch + 1);

    label nInternalFace = nPatchWithFace_ - 1;
    label startFace = nInternalFace;
    label entryi = 0;
    for (const entry& e : patchEntries_)
    {
        // Re-create boundary types, but reset nFaces and startFace settings
        dictionary patchDict = e.dict();
        const word& patchName = e.keyword();

        DebugPout << "Setting " << patchName << endl;

        label nFaces0 = patchDict.get<label>("nFaces");

        if (nFaces0)
        {
            // Only set to 4 faces if there were faces in the original patch
            nFaces0 = 4;
            patchDict.set("nFaces", nFaces0);
        }

        patchDict.set("startFace", startFace);
        patches[entryi] =
            polyPatch::New
            (
                patchName,
                patchDict,
                entryi,
                mesh.boundaryMesh()
            ).ptr();

        ++entryi;
        startFace += nFaces0;
    }

    patches.last() = new emptyPolyPatch
    (
        typeName + ":default",              // name
        2,                                  // number of faces
        nInternalFace + 4*nPatchWithFace_,  // start face
        nPatch - 1,                         // index in boundary list
        mesh.boundaryMesh(),                // polyBoundaryMesh
        emptyPolyPatch::typeName            // patchType
    );

    mesh.addFvPatches(patches);

    if (debug)
    {
        Pout<< "patches:" << nl << mesh.boundaryMesh() << endl;
    }
}


void Foam::simplifiedMeshes::columnFvMeshInfo::initialiseZones(fvMesh& mesh)
{
    if (createFromMesh_)
    {
        // Initialise the zones
        initialiseZone<pointZoneMesh>
        (
            "point",
            localInstance_,
            mesh.pointZones()
        );
        initialiseZone<faceZoneMesh>("face", localInstance_, mesh.faceZones());
        initialiseZone<cellZoneMesh>("cell", localInstance_, mesh.cellZones());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedMeshes::columnFvMeshInfo::columnFvMeshInfo
(
    const Time& runTime,
    const word& regionName
)
:
    regionName_(regionName),
    localInstance_
    (
        runTime.findInstance
        (
            polyMesh::meshDir(regionName_),
            "boundary",
            IOobject::READ_IF_PRESENT
        )
    ),
    createFromMesh_(false),
    points1D_(),
    faces1D_(),
    owner1D_(),
    neighbour1D_(),
    patchEntries_(),
    nPatchWithFace_(0)
{
    initialise(runTime);

    // Dummy zones and sets created on demand
    // Note: zones can be updated post-construction
    cellZoneMesh::disallowGenericZones = 1;
    topoSet::disallowGenericSets = 1;
}


Foam::simplifiedMeshes::columnFvMesh::columnFvMesh
(
    const Time& runTime,
    const word& regionName
)
:
    columnFvMeshInfo(runTime, regionName),
    simplifiedFvMesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime,
            IOobject::NO_READ, // Do not read any existing mesh
            IOobject::NO_WRITE
        ),
        std::move(points1D_),
        std::move(faces1D_),
        std::move(owner1D_),
        std::move(neighbour1D_)
    )
{
    // Workaround to read fvSchemes and fvSolution after setting NO_READ
    // when creating the mesh
    {
        fvSchemes::readOpt(IOobject::MUST_READ);
        fvSchemes::read();
        fvSolution::readOpt(IOobject::MUST_READ);
        fvSolution::read();
    }

    // Add the patches
    addLocalPatches(*this);

    // Add the zones if constructed from mesh
    initialiseZones(*this);

    if (debug)
    {
        setInstance(runTime.timeName());
        objectRegistry::write();
    }
}


// ************************************************************************* //
