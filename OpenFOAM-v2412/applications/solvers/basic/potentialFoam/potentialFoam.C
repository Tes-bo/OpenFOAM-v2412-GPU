/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    potentialFoam

Group
    grpBasicSolvers

Description
    Potential flow solver which solves for the velocity potential, to
    calculate the flux-field, from which the velocity field is obtained by
    reconstructing the flux.

    \heading Solver details
    The potential flow solution is typically employed to generate initial fields
    for full Navier-Stokes codes.  The flow is evolved using the equation:

    \f[
        \laplacian \Phi = \div(\vec{U})
    \f]

    Where:
    \vartable
        \Phi      | Velocity potential [m2/s]
        \vec{U}   | Velocity [m/s]
    \endvartable

    The corresponding pressure field could be calculated from the divergence
    of the Euler equation:

    \f[
        \laplacian p + \div(\div(\vec{U}\otimes\vec{U})) = 0
    \f]

    but this generates excessive pressure variation in regions of large
    velocity gradient normal to the flow direction.  A better option is to
    calculate the pressure field corresponding to velocity variation along the
    stream-lines:

    \f[
        \laplacian p + \div(\vec{F}\cdot\div(\vec{U}\otimes\vec{U})) = 0
    \f]
    where the flow direction tensor \f$\vec{F}\f$ is obtained from
    \f[
        \vec{F} = \hat{\vec{U}}\otimes\hat{\vec{U}}
    \f]

    \heading Required fields
    \plaintable
        U         | Velocity [m/s]
    \endplaintable

    \heading Optional fields
    \plaintable
        p         | Kinematic pressure [m2/s2]
        Phi       | Velocity potential [m2/s]
                  | Generated from p (if present) or U if not present
    \endplaintable

    \heading Options
    \plaintable
        -writep   | write the Euler pressure
        -writephi | Write the final volumetric flux
        -writePhi | Write the final velocity potential
        -initialiseUBCs | Update the velocity boundaries before solving for Phi
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Potential flow solver which solves for the velocity potential"
    );

    argList::addOption
    (
        "pName",
        "pName",
        "Name of the pressure field"
    );

    argList::addBoolOption
    (
        "initialiseUBCs",
        "Initialise U boundary conditions"
    );

    argList::addBoolOption
    (
        "writephi",
        "Write the final volumetric flux field"
    );

    argList::addBoolOption
    (
        "writePhi",
        "Write the final velocity potential field"
    );

    argList::addBoolOption
    (
        "writep",
        "Calculate and write the Euler pressure field"
    );

    argList::addBoolOption
    (
        "withFunctionObjects",
        "Execute functionObjects"
    );

    #include "addRegionOption.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl potentialFlow(mesh, "potentialFlow");

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating potential flow" << endl;

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    MRF.makeRelative(phi);
    adjustPhi(phi, U, p);

    // Non-orthogonal velocity potential corrector loop
    while (potentialFlow.correctNonOrthogonal())
    {
        fvScalarMatrix PhiEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
         ==
            fvc::div(phi)
        );

        PhiEqn.setReference(PhiRefCell, PhiRefValue);
        PhiEqn.solve();

        if (potentialFlow.finalNonOrthogonalIter())
        {
            phi -= PhiEqn.flux();
        }
    }

    MRF.makeAbsolute(phi);

    Info<< "Continuity error = "
        << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
        << endl;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Info<< "Interpolated velocity error = "
        << (sqrt(sum(sqr(fvc::flux(U) - phi)))/sum(mesh.magSf())).value()
        << endl;

    // Write U
    U.write();

    // Optionally write the volumetric flux, phi
    if (args.found("writephi"))
    {
        phi.write();
    }

    // Optionally write velocity potential, Phi
    if (args.found("writePhi"))
    {
        Phi.write();
    }

    // Calculate the pressure field from the Euler equation
    if (args.found("writep"))
    {
        Info<< nl << "Calculating approximate pressure field" << endl;

        label pRefCell = 0;
        scalar pRefValue = 0.0;
        setRefCell
        (
            p,
            potentialFlow.dict(),
            pRefCell,
            pRefValue
        );

        // Calculate the flow-direction filter tensor
        volScalarField magSqrU(magSqr(U));
        volSymmTensorField F(sqr(U)/(magSqrU + SMALL*average(magSqrU)));

        // Calculate the divergence of the flow-direction filtered div(U*U)
        // Filtering with the flow-direction generates a more reasonable
        // pressure distribution in regions of high velocity gradient in the
        // direction of the flow
        volScalarField divDivUU
        (
            fvc::div
            (
                F & fvc::div(phi, U),
                "div(div(phi,U))"
            )
        );

        // Solve a Poisson equation for the approximate pressure
        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(p) + divDivUU
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();
        }

        p.write();
    }

    runTime.functionObjects().end();

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
