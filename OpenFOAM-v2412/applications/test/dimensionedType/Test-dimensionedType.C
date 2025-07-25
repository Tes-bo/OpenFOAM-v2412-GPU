/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "primitiveEntry.H"
#include "dimensionedTypes.H"
#include "SpanStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    dimensionedTensor dt("dt", dimLength, tensor(0, 1, 2, 3, 4, 5, 6, 7, 8));

    Info<< "dt.component(tensor::XX): " << dt.component(tensor::XX) << endl;

    dimensionedScalar ds("ds", dimTime, 1.0);

    Info<< "ds*dt dt*ds: " << ds*dt << " " << dt*ds << endl;

    dimensionedTensor dt2("dt2", dimLength, tensor(1, 1, 2, 3, 4, 5, 6, 7, 8));

    Info<< nl
        << "cmptMultiply(dt, dt2): " << cmptMultiply(dt, dt2) << nl
        << "cmptDivide(dt, dt2): " << cmptDivide(dt, dt2) << endl;

    {
        string dimString("[Pa m^2 s^-2]");
        ICharStream is("[Pa m^2 s^-2]");

        Pout<< nl;
        Pout<< "dimensionSet construct from (is) with contents "
            << dimString << " = " << dimensionSet(is)
            << endl;

        is.rewind();

        dimensionSet dset(dimless);
        is >> dset;
        Pout<< "dimensionSet read:" << dset << endl;
    }

    {
        Pout<< nl
            << "construct from (is) = "
            << dimensionedScalar(ICharStream("bla [Pa mm^2 s^-2] 3.0")())
            << endl;
        Pout<< "construct from (name,is) = "
            <<  dimensionedScalar
                (
                    "ABC",
                    ICharStream("[Pa mm^2 s^-2] 3.0")()
                ) << endl;
        Pout<< "construct from (name,dims,is) = "
            <<  dimensionedScalar
                (
                    "ABC",
                    dimLength,
                    ICharStream("bla [mm] 3.0")()
                ) << endl;
        {
            ICharStream is("bla [mm] 3.0");
            dimensionedScalar ds;
            is >> ds;
            Pout<< "read:" << ds << endl;

            Info<< "writeEntry:" << nl;

            Info<< "abc> "; ds.writeEntry("abc", Info);
            Info<< endl;

            Info<< "bla> "; ds.writeEntry("bla", Info);
            Info<< endl;
        }
    }

    {
        dimensionedLabel a("min", dimTime, -10);
        dimensionedLabel b("max", dimTime, 100);
        scalar t = 0.5;

        Info<< "lerp of" << nl
            << "    " << a << nl
            << "    " << b << nl
            << "  = " << lerp(a, b, t) << nl
            << " vs " << ((1-t)*a + t*b) << nl
            << nl;
    }

    Pout<< "zero scalar (time): " << dimensionedScalar(dimTime) << endl;
    Pout<< "zero vector: " << dimensionedVector(dimLength) << endl;
    Pout<< "zero tensor: " << dimensionedTensor(dimLength) << endl;


    dictionary dict;
    {
        dict.add("test1", scalar(10));
        dict.add("test2a", scalar(21));
        dict.add("test5", dimensionedScalar("", 50));
        dict.add("carp1", dimensionedScalar("test1", 11));
        // This will fail to tokenize:
        // dict.add("test5", dimensionedScalar(50));
    }

    Info<< nl << "Get from dictionary: " << dict << nl;

    Info<< "test1 : " << dimensionedScalar("test1", dict) << nl;
    Info<< "test2 : " << dimensionedScalar("test2", dimless, 20, dict) << nl;
    Info<< "test2a : " << dimensionedScalar("test2a", dimless, 20, dict) << nl;
    Info<< "test3 : "
        << dimensionedScalar::getOrDefault("test3", dict, 30) << nl;

    Info<< "test4 : "
        << dimensionedScalar::getOrAddToDict("test4", dict, 40) << nl;

    Info<< "test5 : "
        << dimensionedScalar::getOrAddToDict("test5", dict, -50) << nl;

    // Deprecated
    Info<< "Deprecated constructors" << nl;
    Info<< "carp : "
        << dimensionedScalar(dict.lookup("carp1")) << nl;

    Info<< "carp : "
        << dimensionedScalar("other", dict.lookup("test5")) << nl;

    Info<< "carp : "
        << dimensionedScalar("carp", dimless, dict.lookup("carp1")) << nl;

    Info<< "alt : "
        << dimensionedScalar("myName", dimless, dict, "carp1") << nl;

    Info<< "alt : "
        << dimensionedScalar("myName", dimless, dict, "test5") << nl;

    {
        dimensionedScalar scalar1("myName", dimless, Zero);

        scalar1.read("test5", dict);
        Info<< "read in : " << scalar1 << nl;

        scalar1.readIfPresent("test4", dict);
        Info<< "read in : " << scalar1 << nl;

        scalar1.readIfPresent("test5", dict);
        Info<< "read in : " << scalar1 << nl;
    }

    Info<< nl << "Dictionary is now: " << dict << nl;


    {
        primitiveEntry entry1("scalar1", token(15.0));

        // The same:
        // primitiveEntry entry1("scalar1", ITstream::parse(" 15.0 "));

        // This fails (as it should):
        // primitiveEntry entry1("scalar1", ITstream::parse(" 15.0 25.0 "));

        dimensionedScalar ds1(entry1);

        Info<< "construct from entry: "
            << entry1 << nl
            << "    = " << ds1 << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
