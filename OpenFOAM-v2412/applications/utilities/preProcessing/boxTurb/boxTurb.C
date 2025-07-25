/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    boxTurb

Group
    grpPreProcessingUtilities

Description
    Create a box of divergence-free turbulence conforming to a given
    energy spectrum.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "graph.H"
#include "OFstream.H"
#include "Kmesh.H"
#include "turbGen.H"
#include "calcEk.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a box of divergence-free turbulence conforming to a given"
        " energy spectrum"
    );

    argList::noParallel();
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFields.H"
    #include "readBoxTurbDict.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Kmesh K(mesh);

    turbGen Ugen(K, Ea, k0);

    U.primitiveFieldRef() = Ugen.U();
    U.correctBoundaryConditions();

    Info<< "k("
         << runTime.timeName()
         << ") = "
         << 3.0/2.0*average(magSqr(U)).value() << endl;

    U.write();

    calcEk(U, K).write
    (
        runTime.path()/"graphs"/runTime.timeName(),
        "Ek",
        runTime.graphFormat()
    );

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
