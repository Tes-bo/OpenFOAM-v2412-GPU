/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2016 OpenCFD Ltd.
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
    XiFoam

Group
    grpCombustionSolvers grpMovingMeshSolvers

Description
    Solver for compressible premixed/partially-premixed combustion with
    turbulence modelling.

    Combusting RANS code using the b-Xi two-equation model.
    Xi may be obtained by either the solution of the Xi transport
    equation or from an algebraic expression.  Both approaches are
    based on Gulder's flame speed correlation which has been shown
    to be appropriate by comparison with the results from the
    spectral model.

    Strain effects are encorporated directly into the Xi equation
    but not in the algebraic approximation.  Further work need to be
    done on this issue, particularly regarding the enhanced removal rate
    caused by flame compression.  Analysis using results of the spectral
    model will be required.

    For cases involving very lean Propane flames or other flames which are
    very strain-sensitive, a transport equation for the laminar flame
    speed is present.  This equation is derived using heuristic arguments
    involving the strain time scale and the strain-rate at extinction.
    the transport velocity is the same as that for the Xi equation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiuReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "laminarFlameSpeed.H"
#include "ignition.H"
#include "Switch.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for compressible premixed/partially-premixed combustion with"
        " turbulence modelling."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "readCombustionProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createRhoUf.H"
    #include "createControls.H"
    #include "initContinuityErrs.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "createTimeControls.H"

        {
            // Store divrhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
            volScalarField divrhoU
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );

            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"

            ++runTime;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Store momentum to set rhoUf for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            mesh.update();

            if (mesh.changing() && correctPhi)
            {
                // Calculate absolute flux from the mapped surface velocity
                phi = mesh.Sf() & rhoUf;

                #include "correctPhi.H"

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        #include "rhoEqn.H"
        Info<< "rho min/max : " << min(rho).value() << " " << max(rho).value()
            << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            #include "ftEqn.H"
            #include "bEqn.H"
            #include "EauEqn.H"
            #include "EaEqn.H"

            if (!ign.ignited())
            {
                thermo.heu() == thermo.he();
            }

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

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
