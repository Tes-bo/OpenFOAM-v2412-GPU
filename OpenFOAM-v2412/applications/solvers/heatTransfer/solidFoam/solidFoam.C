/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    solidFoam

Group
    grpHeatTransferSolvers

Description
    Solver for energy transport and thermodynamics on a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "dummyCourantNo.H"
#include "solidRegionDiffNo.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for energy transport and thermodynamics on a solid"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEvolving thermodynamics\n" << endl;

    if (mesh.solutionDict().found("SIMPLE"))
    {
        simpleControl simple(mesh);

        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            while (simple.correctNonOrthogonal())
            {
                #include "hEqn.H"
            }

            runTime.write();

            runTime.printExecutionTime(Info);
        }
    }
    else
    {
        pimpleControl pimple(mesh);

        #include "createDyMControls.H"

        while (runTime.run())
        {
            #include "readDyMControls.H"
            #include "readSolidTimeControls.H"

            #include "solidDiffusionNo.H"
            #include "setMultiRegionDeltaT.H"

            ++runTime;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    mesh.controlledUpdate();

                    if (mesh.changing() && checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                while (pimple.correct())
                {
                    #include "hEqn.H"
                }
            }

            runTime.write();
            runTime.printExecutionTime(Info);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
