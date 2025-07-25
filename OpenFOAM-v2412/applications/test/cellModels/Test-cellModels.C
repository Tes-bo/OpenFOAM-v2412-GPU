/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Test-cellModels

Description
    Print information about known cellModels

\*---------------------------------------------------------------------------*/

#include "cellModel.H"
#include "cellModeller.H"

using namespace Foam;

void printInfo(const cellModel* mdl)
{
    if (mdl)
    {
        Info<< *mdl << endl;
    }
    else
    {
        Info<< "nullptr" << endl;
    }
}


void printInfo(const cellModel::modelType type)
{
    Info<< cellModel::modelNames[type] << " = ";
    printInfo(cellModel::ptr(type));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<<"lookup by enum" << nl
        <<"=========================" << endl;

    printInfo(cellModel::UNKNOWN);
    printInfo(cellModel::HEX);
    printInfo(cellModel::WEDGE);
    printInfo(cellModel::PRISM);
    printInfo(cellModel::PYR);
    printInfo(cellModel::TET);
    printInfo(cellModel::SPLITHEX);
    printInfo(cellModel::TETWEDGE);


    Info<<"lookup by name" << nl
        <<"=========================" << endl;

    printInfo(cellModel::ptr("tet"));

    Info<<"lookup by index" << nl
        <<"=========================" << endl;

    printInfo(cellModel::ptr(7));

    // Compatibility mode
    Info<<"cellModeller::lookup (compatibility)" << nl
        <<"=========================" << endl;

    printInfo(cellModeller::lookup("tet"));

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
