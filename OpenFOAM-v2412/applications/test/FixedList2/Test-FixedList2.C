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
    Test-FixedList2

Description
    Test speeds, usability of some List/FixedList operations

See also
    Foam::FixedList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "FixedList.H"
#include "labelList.H"
#include "vectorList.H"
#include "ListOps.H"
#include "IFstream.H"
#include "OFstream.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

template<class ListType>
void runSwapTest
(
    const label nLoops,
    ListType& list1,
    ListType& list2
)
{
    cpuTime timer;

    Info<<"Swapping fixed lists with " << list1.size() << " elements\n";

    Info<< "input 1: " << list1.first() << nl;
    Info<< "input 2: " << list2.first() << nl;

    // Should be zero, since this is a compile-time value

    Info<< "Perform " << nLoops << " swaps..." << nl;

    for (label iLoop = 0; iLoop < nLoops; ++iLoop)
    {
        Foam::Swap(list1, list2);
    }

    Info<< "output 1: " << list1.first() << nl;
    Info<< "output 2: " << list2.first() << nl;

    Info<< "Operation took"
        << "  " << timer.cpuTimeIncrement() << " s\n\n";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("label");
    argList::addBoolOption("float");
    argList::addBoolOption("vector");
    argList::addBoolOption("labelList");
    argList::addBoolOption("vectorList");
    argList::addBoolOption("fixedLabel");
    argList::addBoolOption("fixedLabelList");

    argList args(argc, argv);

    if (args.options().empty())
    {
        Info<< nl << "Specify an option! " << nl << endl;
    }

    if (args.found("label"))
    {
        FixedList<label, 100000> list1(1);
        FixedList<label, 100000> list2(0);

        runSwapTest(1000001, list1, list2);
    }

    if (args.found("float"))
    {
        FixedList<double, 100000> list1(1.0);
        FixedList<double, 100000> list2(0.0);

        runSwapTest(1000001, list1, list2);
    }

    if (args.found("vector"))
    {
        FixedList<vector, 100000> list1(vector::one);
        FixedList<vector, 100000> list2(vector::zero);

        runSwapTest(100001, list1, list2);
    }

    if (args.found("labelList"))
    {
        typedef labelList testType;
        testType initVal(500);

        initVal = 0;
        FixedList<testType, 1000> list1(initVal);

        initVal = 1;
        FixedList<testType, 1000> list2(initVal);

        runSwapTest(100001, list1, list2);
    }

    if (args.found("vectorList"))
    {
        typedef vectorList testType;
        testType initVal(500);

        initVal = vector::zero;
        FixedList<testType, 1000> list1(initVal);

        initVal = vector::one;
        FixedList<testType, 1000> list2(initVal);

        runSwapTest(100001, list1, list2);
    }

    if (args.found("fixedLabel"))
    {
        typedef FixedList<label,1000> testType;

        testType initVal;

        initVal = 0;
        FixedList<testType, 1000> list1(initVal);

        initVal = 1;
        FixedList<testType, 1000> list2(initVal);

        runSwapTest(100001, list1, list2);
    }

    if (args.found("fixedLabelList"))
    {
        typedef labelList testType;
        typedef FixedList<testType,10> containerType;

        testType tinitVal(500);
        containerType initVal;

        tinitVal = 0;
        initVal = tinitVal;
        FixedList<containerType, 1000> list1(initVal);

        tinitVal = 1;
        initVal = tinitVal;
        FixedList<containerType, 1000> list2(initVal);

        runSwapTest(10001, list1, list2);
    }

    Info<< nl << "Done" << nl << endl;
    return 0;
}


// ************************************************************************* //
