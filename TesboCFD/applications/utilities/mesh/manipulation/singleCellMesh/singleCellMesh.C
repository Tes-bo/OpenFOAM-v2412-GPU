/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    singleCellMesh

Group
    grpMeshManipulationUtilities

Description
    Reads all fields and maps them to a mesh with all internal faces removed
    (singleCellFvMesh) which gets written to region "singleCell".

    Used to generate mesh and fields that can be used for boundary-only data.
    Might easily result in illegal mesh though so only look at boundaries
    in paraview.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "ReadFields.H"
#include "singleCellFvMesh.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Name of region to create
const word singleCellName = "singleCell";


template<class GeoField>
wordList interpolateFields
(
    const singleCellFvMesh& subsetter,
    const PtrList<GeoField>& flds
)
{
    wordList names(flds.size());

    forAll(flds, i)
    {
        GeoField* subFld = subsetter.interpolate(flds[i]).ptr();

        subFld->writeOpt(IOobject::AUTO_WRITE);
        subFld->store();

        names[i] = subFld->name();
    }

    return names;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Map fields to a mesh with all internal faces removed"
        " (singleCellFvMesh) which gets written to region 'singleCell'"
    );

    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    if (regionName == singleCellName)
    {
        FatalErrorInFunction
            << "Cannot convert region " << regionName
            << " since result would overwrite it. Please rename your region."
            << exit(FatalError);
    }

    // Create the mesh
    Info<< "Creating singleCell mesh" << nl << endl;
    autoPtr<singleCellFvMesh> scMesh
    (
        new singleCellFvMesh
        (
            IOobject
            (
                singleCellName,
                mesh.pointsInstance(),
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh
        )
    );
    scMesh().setInstance(mesh.pointsInstance());

    // For convenience create any fv* files
    if (!exists(scMesh().fvSolution::objectPath()))
    {
        mkDir(scMesh().fvSolution::path());
        ln("../fvSolution", scMesh().fvSolution::objectPath());
    }
    if (!exists(scMesh().fvSchemes::objectPath()))
    {
        mkDir(scMesh().fvSolution::path());
        ln("../fvSchemes", scMesh().fvSchemes::objectPath());
    }


    // List of stored objects to clear prior to reading
    DynamicList<word> storedObjects;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< nl << "Time = " << runTime.timeName() << endl;

        // Purge any previously interpolated fields
        if (!storedObjects.empty())
        {
            static_cast<objectRegistry&>(scMesh()).erase(storedObjects);
            storedObjects.clear();
        }

        // Check for new mesh
        if (mesh.readUpdate() != polyMesh::UNCHANGED)
        {
            Info<< "Detected changed mesh. Recreating singleCell mesh." << endl;
            scMesh.reset(nullptr);  // first free any stored objects
            scMesh.reset
            (
                new singleCellFvMesh
                (
                    IOobject
                    (
                        singleCellName,
                        mesh.pointsInstance(),
                        runTime,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    mesh
                )
            );
        }

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());
        storedObjects.reserve(objects.size());

        // Read fvMesh fields, map and store the interpolated fields
        #define doLocalCode(GeoField)                                         \
        {                                                                     \
            PtrList<GeoField> flds;                                           \
            ReadFields(mesh, objects, flds);                                  \
            storedObjects.push_back(interpolateFields(scMesh(), flds));       \
        }

        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);

        #undef doLocalCode

        // Write
        Info<< "Writing " << singleCellName
            << " mesh/fields to time " << runTime.timeName() << endl;
        scMesh().write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
