/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "predicates.H"
#include "dynamicFvMesh.H"
#include "dynamicOversetFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createNamedDynamicFvMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCorrecting boundary conditions on " << T.name() << nl << endl;

    ++runTime;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    Info<< "Reading : ";
    runTime.printExecutionTime(Info);

    mesh.update();

    Info<< "Overset calculation : ";
    runTime.printExecutionTime(Info);

    if (false)
    {
        // Test correctBoundaryConditions

        // Change the internalField
        component(T.ref(), mesh.C(), 0);
        component(T.ref(), mesh.C(), 1);
        // Interpolate + halo swap
        T.correctBoundaryConditions();
        // T.boundaryFieldRef().evaluate_if(predicates::always{});
        // Check halo swap
        dynamicOversetFvMesh::checkCoupledBC(T);
    }
    if (true)
    {
        // Test solving
        fvScalarMatrix TEqn(fvm::laplacian(DT, T));
        TEqn.solve();
    }

    runTime.write();

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
