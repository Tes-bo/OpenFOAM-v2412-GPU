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
    Test-vtmWriter

Description
    Basic functionality tests for vtk::vtmWriter

\*---------------------------------------------------------------------------*/

#include "foamVtmWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    vtk::vtmWriter writer1;
    {
        fileName base = "region1_0001";

        writer1.beginBlock("internal");
        writer1.append_vtu
        (
            base/"internal"
        );
        writer1.endBlock("internal");

        {
            writer1.beginBlock("boundary");
            writer1.append_vtp
            (
                base/"patch0"
            );
            writer1.append("");  // bad entry
            writer1.append_vtp
            (
                base/"patch1"
            );
            writer1.append_vtp
            (
                base/"patch2"
            );
        }

        writer1.endBlock("boundary");

        {
            writer1.beginBlock("empty");
            writer1.endBlock("empty");
        }
        {
            writer1.beginBlock("dangling1");
            writer1.beginBlock("dangling2");
        }
    }

    Info<< nl << "vtm information" << nl;
    writer1.dump(Info),
    Info<< nl;

//    writer1.repair();
//
//   Info<< nl << "vtm information - after repair" << nl;
//    writer1.dump(Info),
//    Info<< nl;

    writer1.repair(true);

//    Info<< nl << "vtm information - after repair(collapse)" << nl;
//    writer1.dump(Info),
//    Info<< nl;
//
//    Info<< nl << "vtm information - after repair(collapse)" << nl;
//    writer1.dump(Info),
//    Info<< nl;

    Info<< nl << "Write to file" << nl;
    writer1.write("vtmWriter1.vtm");


    vtk::vtmWriter writer2;
    {
        fileName base = "region2_0001";

        writer2.beginBlock("internal");
        writer2.append_vtu
        (
            base/"internal"
        );
        writer2.endBlock("internal");

        {
            writer2.beginBlock("boundary");
            writer2.append_vtp
            (
                base/"patch0"
            );
            writer2.append("");  // bad entry
            writer2.append_vtp
            (
                base/"patch1"
            );
            writer2.append_vtp
            (
                base/"patch2"
            );
        }

        writer2.endBlock("boundary");

        // These should be automatically skipped
        writer2.endBlock();
        writer2.endBlock();
        writer2.endBlock();
        writer2.endBlock();
    }

    Info<< nl << "vtm information" << nl;
    writer2.dump(Info);

    writer2.repair(true);


    vtk::vtmWriter writer3;

    writer3.add("some-region1", writer1);
    writer3.add("some-region2", writer2);

    Info<< nl << "Combined:" << nl;
    writer3.dump(Info);

    return 0;
}


// ************************************************************************* //
