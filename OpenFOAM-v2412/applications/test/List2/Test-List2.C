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
    Test-List2

Description
    Test speeds, usability of some List/FixedList operations

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
void runResizeTest
(
    const label nLoops,
    ListType& list,
    std::initializer_list<label> sizes
)
{
    cpuTime timer;

    const label size0 = list.size();
    const auto val = list.first();

    Info<<"Resize list(" << list.size() << ") to";

    for (auto len : sizes)
    {
        Info<< " " << len;
    }
    Info<< nl;

    Info<< "Perform " << nLoops << " times..." << nl;
    for (label iLoop = 0; iLoop < nLoops; ++iLoop)
    {
        list.setSize(size0, val);

        for (auto len : sizes)
        {
            list.setSize(len, val);
        }
    }

    Info<< "Operation took"
        << "  " << timer.cpuTimeIncrement() << " s\n\n";
}


template<class ListType>
void runOrderingTest(const label nLoops, const ListType& list)
{
    if (true)
    {
        cpuTime timer;

        float total = 0;

        Info<<"forAll - perform " << nLoops << " times..." << nl;
        for (label iLoop = 0; iLoop < nLoops; ++iLoop)
        {
            float sum = 0;
            forAll(list, i)
            {
                sum += list[i];
            }

            total += sum;
        }

        Info<< "Operation (sum " << total << ") took"
            << "  " << timer.cpuTimeIncrement() << " s\n\n";
    }

    if (true)
    {
        cpuTime timer;

        float total = 0;

        Info<<"reverse pointer loop - perform " << nLoops << " times..." << nl;
        for (label iLoop = 0; iLoop < nLoops; ++iLoop)
        {
            float sum = 0;

            const typename ListType::value_type* __restrict__ fp
                = (list).end();

            label i = (list).size();
            while (i--)
            {
                sum += (*--fp);
            }

            total += sum;
        }

        Info<< "Operation (sum " << total << ") took"
            << "  " << timer.cpuTimeIncrement() << " s\n\n";
    }

    if (true)
    {
        cpuTime timer;

        float total = 0;

        Info<<"forward pointer loop - perform " << nLoops << " times..." << nl;
        for (label iLoop = 0; iLoop < nLoops; ++iLoop)
        {
            float sum = 0;

            const typename ListType::value_type* __restrict__ fp
                = (list).begin();

            label i = (list).size();
            while (i--)
            {
                sum += (*fp++);
            }

            total += sum;
        }

        Info<< "Operation (sum " << total << ") took"
            << "  " << timer.cpuTimeIncrement() << " s\n\n";
    }


    if (true)
    {
        cpuTime timer;

        float total = 0;

        Info<<"for loop - perform " << nLoops << " times..." << nl;
        for (label iLoop = 0; iLoop < nLoops; ++iLoop)
        {
            float sum = 0;

            const typename ListType::value_type* __restrict__ fp
                = (list).begin();

            const label sz = (list).size();
            for (label i=0; i<sz; ++i)
            {
                sum += fp[i];
            }

            total += sum;
        }

        Info<< "Operation (sum " << total << ") took"
            << "  " << timer.cpuTimeIncrement() << " s\n\n";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("label");
    argList::addBoolOption("float");
    argList::addBoolOption("vector");
    argList::addBoolOption("order");
    argList::addBoolOption("labelList");
    argList::addBoolOption("vectorList");

    argList args(argc, argv);

    if (args.options().empty())
    {
        Info<< nl << "Specify an option! " << nl << endl;
    }


    std::initializer_list<label> increments
        = {10000, 20000, 40000, 80000, 160000};

    if (args.found("label"))
    {
        List<label> list(10, 1);

        runResizeTest(100000, list, increments);
    }

    if (args.found("float"))
    {
        List<double> list(10, 1.0);

        runResizeTest(10000, list, increments);
    }

    if (args.found("vector"))
    {
        List<vector> list(10, vector::one);

        runResizeTest(10000, list, increments);
    }

    if (args.found("labelList"))
    {
        typedef labelList testType;
        testType initVal(500, label(1));

        List<testType> list(10, initVal);

        runResizeTest(200, list, increments);
    }

    if (args.found("vectorList"))
    {
        typedef vectorList testType;
        testType initVal(500, vector::one);

        List<testType> list(10, initVal);

        runResizeTest(100, list, increments);
    }

    if (args.found("order"))
    {
        List<label> list(100000000, 1);

        runOrderingTest(100, list);
    }


    Info<< nl << "Done" << nl << endl;
    return 0;
}


// ************************************************************************* //
