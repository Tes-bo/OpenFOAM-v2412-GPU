/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
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
    chtMultiRegionFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow and solid heat"
        " conduction with conjugate heat transfer"
        " between solid and fluid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createCoupledRegions.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }

            if (coupled)
            {
                Info<< "\nSolving energy coupled regions " << endl;
                fvMatrixAssemblyPtr->solve();
                #include "correctThermos.H"

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    if (!frozenFlow)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        // --- PISO loop
                        for (int corr=0; corr<nCorr; corr++)
                        {
                            #include "pEqn.H"
                        }
                        turbulence.correct();
                    }

                    rho = thermo.rho();
                    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                        << max(thermo.T()).value() << endl;
                }

                fvMatrixAssemblyPtr->clear();
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        fvMesh& mesh = fluidRegions[i];

                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "readFluidMultiRegionPIMPLEControls.H"
                        #include "setRegionFluidFields.H"
                        frozenFlow = true;
                        #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        fvMesh& mesh = solidRegions[i];

                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "setRegionSolidFields.H"
                        #include "solveSolid.H"
                    }

                    if (coupled)
                    {
                        Info<< "\nSolving energy coupled regions " << endl;
                        fvMatrixAssemblyPtr->solve();
                        #include "correctThermos.H"

                        forAll(fluidRegions, i)
                        {
                            #include "setRegionFluidFields.H"
                            rho = thermo.rho();
                        }

                        fvMatrixAssemblyPtr->clear();
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
