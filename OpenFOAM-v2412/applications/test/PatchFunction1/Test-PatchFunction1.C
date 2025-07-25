/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test-Function1

Description
    Tests Function1

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PatchFunction1.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary function1Properties
    (
        IOobject
        (
            "function1Properties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    autoPtr<PatchFunction1<vector>> function1
    (
        PatchFunction1<vector>::New
        (
            mesh.boundaryMesh()[0],
            "function1",
            function1Properties
        )
    );

    scalar x0 = function1Properties.get<scalar>("x0");
    scalar x1 = function1Properties.get<scalar>("x1");

    Info<< "Data entry type: " << function1().type() << nl << endl;

    Info<< "Inputs" << nl
        << "    x0 = " << x0 << nl
        << "    x1 = " << x1 << nl
        << endl;

    Info<< "Interpolation" << nl
        << "    f(x0) = " << function1().value(x0) << nl
        << "    f(x1) = " << function1().value(x1) << nl
        << endl;

    Info<< "Integration" << nl
        << "    int(f(x)) lim(x0->x1) = " << function1().integrate(x0, x1) << nl
        << endl;

    volVectorField fld
    (
        IOobject
        (
            "value",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedVector(Zero)
    );

    fld.boundaryFieldRef()[0] == function1().value(x0);
    fld.write();

    return 0;
}


// ************************************************************************* //
