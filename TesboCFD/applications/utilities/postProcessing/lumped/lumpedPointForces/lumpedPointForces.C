/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    lumpedPointForces

Description
    Extract force/moment information from simulation results that
    use the lumped points movement description.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "foamVtkSeriesWriter.H"
#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoFieldType>
autoPtr<GeoFieldType> loadField
(
    const volMesh::Mesh& mesh,
    const IOobject* io
)
{
    if (io && io->isHeaderClass<GeoFieldType>())
    {
        Info<< "Reading " << io->headerClassName()
            << ' ' << io->name() << endl;

        return autoPtr<GeoFieldType>::New
        (
            IOobject
            (
                io->name(),
                io->instance(),
                io->local(),
                io->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                io->registerObject()
            ),
            mesh
        );
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract force/moment information from simulation results that"
        " use the lumped points movement description."
    );

    argList::addBoolOption
    (
        "vtk",
        "Create visualization files of the forces"
    );

    timeSelector::addOptions(true, false);
    argList::noFunctionObjects();  // Never use function objects

    #include "addRegionOption.H"
    #include "setRootCase.H"

    const bool withVTK = args.found("vtk");

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New(mesh);

    if (!movement)
    {
        Info<< "No valid movement found" << endl;
        return 1;
    }

    const label nPatches = lumpedPointTools::setPatchControls(mesh);
    if (!nPatches)
    {
        Info<< "No point patches with lumped movement found" << endl;
        return 2;
    }

    Info<<"Lumped point patch controls set on " << nPatches
        << " patches" << nl;


    vtk::seriesWriter forceSeries;
    List<vector> forces, moments;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate())
        {
            Info<< "    Read new mesh" << nl;
        }

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        // Pressure field
        autoPtr<volScalarField> pField
            = loadField<volScalarField>(mesh, objects.findObject("p"));

        // The forces per zone
        if (movement().forcesAndMoments(mesh, forces, moments))
        {
            Info<<"forces per zone: " << forces << endl;
            Info<<"moments per zone: " << moments << endl;

            if (withVTK && Pstream::master())
            {
                const word outputName =
                    word::printf("forces_%06d.vtp", runTime.timeIndex());

                Info<<"    " << outputName << endl;

                movement().writeForcesAndMomentsVTP
                (
                    outputName,
                    forces,
                    moments
                );

                forceSeries.append(runTime.timeIndex(), outputName);
            }
        }
    }


    // Create file series

    if (forceSeries.size())
    {
        forceSeries.write("forces.vtp");
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
