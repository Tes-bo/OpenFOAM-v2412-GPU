/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test-PrecisionAdaptor

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "primitiveFields.H"
#include "PrecisionAdaptor.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Field<double> content1(8);
    Field<double> content2(8);

    forAll(content1, i)
    {
        content1[i] = 10 * i;
        content2[i] = 10 * i;
    }

    Foam::reverse(content2);

    ConstPrecisionAdaptor<float, double, Field> cadaptor1;

    cadaptor1.set(content1);
    cadaptor1.commit();  // This is a no-op
    Info<< "wrapped: " << cadaptor1() << nl;

    cadaptor1.set(content2);
    Info<< "wrapped: " << cadaptor1() << nl;

    Info<< nl;

    PrecisionAdaptor<float, double, Field> adaptor2;

    adaptor2.set(content1);
    adaptor2.ref() *= 2;
    adaptor2.commit();  // Propagate changes back to input now

    Info<< "modified wrapped: " << adaptor2() << nl;

    adaptor2.set(content2);
    adaptor2.ref() *= 2;
    adaptor2.commit();  // Propagate changes back to input now

    Info<< "modified wrapped: " << adaptor2() << nl;
    Info<< "source: " << content1 << nl;
    Info<< "source: " << content2 << nl;


    content2 *= 2;
    Info<< nl
        << "set with " << content2 << nl;
    Info<< "wrapped was " << adaptor2() << nl;
    adaptor2.set(content2);
    Info<< "wrapped now " << adaptor2() << nl;
    Info<< "source: " << content2 << nl;

    // Can even do this
    Foam::reverse(adaptor2.ref());

    adaptor2.ref() *= 2;
    adaptor2.set(content1);  // implicit commit
    Info<< "updated: " << content2 << nl;

    Info<< nl
        << "input: " << content1 << nl;

    adaptor2.ref() *= 2;
    adaptor2.clear();  // discard
    adaptor2.commit();

    Info<< "unchanged: " << content1 << nl;

    Info<< nl << "Done" << nl << endl;
    return 0;
}


// ************************************************************************* //
