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
    Test-constantields

Description
    Simple compilation tests for constant fields

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "geometricOneField.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    {
        geometricZeroField fld0;
        geometricOneField fld1;

        Info<< "dims 0: " << fld0.dimensions() << nl;
        Info<< "internal 0: " << scalar(fld0.internalField()[0]) << nl;
        Info<< "boundary 0: " << scalar(fld0.boundaryField()[0][0]) << nl;

        Info<< "dims 1: " << fld1.dimensions() << nl;
        Info<< "internal 1: " << scalar(fld1.internalField()[0]) << nl;
        Info<< "boundary 1: " << scalar(fld1.boundaryField()[0][0]) << nl;
    }

    Info<< "\nDone\n" << endl;

    return 0;
}


// ************************************************************************* //
