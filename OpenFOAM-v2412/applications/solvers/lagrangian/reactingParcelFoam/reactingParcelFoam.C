/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2020 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    reactingParcelFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow with a reacting,
    multiphase particle cloud, and surface film modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "surfaceFilmModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cloudMacros.H"

#ifndef CLOUD_BASE_TYPE
    #define CLOUD_BASE_TYPE ReactingMultiphase
    #define CLOUD_BASE_TYPE_NAME "reacting"
#endif

#include CLOUD_INCLUDE_FILE(CLOUD_BASE_TYPE)
#define basicReactingTypeCloud CLOUD_TYPE(CLOUD_BASE_TYPE)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for compressible, turbulent flow"
        " with reacting, multiphase particle clouds"
        " and surface film modelling."
    );

    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRegionControls.H"
    #include "initContinuityErrs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh
        // so that it can be mapped and used in correctPhi
        // to ensure the corrected phi has the same divergence
        autoPtr<volScalarField> divrhoU;
        if (solvePrimaryRegion && correctPhi)
        {
            divrhoU.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Store momentum to set rhoUf for introduced faces.
        autoPtr<volVectorField> rhoU;
        if (solvePrimaryRegion && rhoUf.valid())
        {
            rhoU.reset(new volVectorField("rhoU", rho*U));
        }

        // Store the particle positions
        parcels.storeGlobalPositions();

        // Do any mesh changes
        mesh.update();

        if (solvePrimaryRegion && mesh.changing())
        {
            gh = (g & mesh.C()) - ghRef;
            ghf = (g & mesh.Cf()) - ghRef;

            MRF.update();

            if (correctPhi)
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phi = mesh.Sf() & rhoUf();

                #include "../../compressible/rhoPimpleFoam/correctPhi.H"

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        parcels.evolve();
        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
            if (pimple.nCorrPIMPLE() <= 1)
            {
                #include "rhoEqn.H"
            }

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
