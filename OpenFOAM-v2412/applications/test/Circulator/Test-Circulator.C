/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Test-Circulator

Description
    Tests for Circulator and ConstCirculator

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "ListOps.H"
#include "face.H"
#include "Circulator.H"


using namespace Foam;

void printCompare(const face& a, const face& b)
{
    Info<< "Compare " << a << " and " << b
        << " Match = " << face::compare(a, b) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test the implementation of a circular iterator" << nl << endl;

    Info<< "Test const circulator. First go forwards, then backwards."
        << nl << endl;

    face f(identity(4));

    // ConstCirculator<face> foo;
    // Info<< "size: " << foo.size() << nl;

    ConstCirculator<face> cStart(f);

    if (!cStart.empty())
    {
        do
        {
            Info<< "Iterate forwards over face (prev/curr/next) : "
                << cStart.prev() << " / "
                << cStart.curr() << " / "
                << cStart.next()
            << endl;
        } while (cStart.circulate(CirculatorBase::CLOCKWISE));
    }

    if (!cStart.empty())
    {
        do
        {
            Info<< "Iterate backwards over face : " << cStart() << endl;

        } while (cStart.circulate(CirculatorBase::ANTICLOCKWISE));
    }


    Info<< nl << nl << "Test non-const circulator" << nl << endl;

    Circulator<face> cStart2(f);

    Info<< "Face before : " << f << endl;

    if (cStart2.size()) do
    {
        Info<< "Iterate forwards over face (prev/curr/next) : "
            << cStart2.prev() << " / " << cStart2() << " / " << cStart2.next()
            << endl;

    } while (cStart2.circulate(CirculatorBase::CLOCKWISE));

    if (!cStart2.empty())
    {
        do
        {
            Info<< "Iterate forwards over face, adding 1 to each element : "
                << cStart2();

            cStart2() += 1;

            Info<< " -> " << cStart2() << endl;
        } while (cStart2.circulate(CirculatorBase::CLOCKWISE));
    }

    Info<< "Face after : " << f << endl;


    Info<< nl << nl << "Compare two faces: " << endl;
    face a(identity(5));
    printCompare(a, a);

    face b(reverseList(a));
    printCompare(a, b);

    face c(a);
    c[4] = 3;
    printCompare(a, c);

    face d(rotateList(a, 2));
    printCompare(a, d);

    face g(labelList(5, 1));
    face h(g);
    printCompare(g, h);

    g[0] = 2;
    h[3] = 2;
    printCompare(g, h);

    g[4] = 3;
    h[4] = 3;
    printCompare(g, h);

    face face1(identity(1));
    printCompare(face1, face1);

    face face2(identity(1, 2));
    printCompare(face1, face2);

    Info<< nl << nl << "Zero face" << nl << endl;

    face fZero;
    Circulator<face> cZero(fZero);

    if (!cZero.empty())
    {
        do
        {
            Info<< "Iterate forwards over face : " << cZero() << endl;

        } while (cZero.circulate(CirculatorBase::CLOCKWISE));
    }

    fZero = face(identity(5));

    // circulator was invalidated so reset
    cZero = Circulator<face>(fZero);

    do
    {
        Info<< "Iterate forwards over face : " << cZero() << endl;

    } while (cZero.circulate(CirculatorBase::CLOCKWISE));


    Info<< nl << nl << "Simultaneously go forwards/backwards over face " << f
        << nl << endl;

    ConstCirculator<face> circForward(f);
    ConstCirculator<face> circBackward(f);

    if (circForward.size() && circBackward.size())
    {
        do
        {
            Info<< "Iterate over face forwards : " << circForward()
                << ", backwards : " << circBackward() << endl;
        }
        while
        (
            circForward.circulate(CirculatorBase::CLOCKWISE),
            circBackward.circulate(CirculatorBase::ANTICLOCKWISE)
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
