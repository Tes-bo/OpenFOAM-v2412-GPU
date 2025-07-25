/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Description
    Test label ranges
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelRanges.H"

using namespace Foam;

void printInfo(const labelRange& range)
{
    Info<< "min     " << range.min() << nl
        << "max     " << range.max() << nl
        << "end     " << range.end_value() << nl
        << "begin/end " << *range.cbegin() << ' ' << *range.cend() << nl;

    // Info<< "rbegin rend " << *range.rbegin() << ' ' << *range.rend() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addArgument("start size .. startN sizeN");
    argList::addVerboseOption("enable labelRange::debug");
    argList::addNote
    (
        "The default is to add ranges, use 'add' and 'del' to toggle\n\n"
        "Eg, 0 10 30 10 del 20 15"
    );

    argList args(argc, argv, false, true);

    if (args.verbose())
    {
        labelRange::debug = 1;
    }

    Info<< nl;
    {
        labelRange range(5, 10);
        Info<< "identity: " << identity(range) << nl;
    }

    {
        Info<< "test sorting" << endl;
        labelRanges list1(10);

        list1.emplace_back(25, 8);
        list1.emplace_back(8);
        list1.emplace_back(15, 5);
        list1.emplace_back(50, -10, true);

        // Move construct
        labelRanges ranges(std::move(list1));
        if (!list1.empty())
        {
            Info<< "Move construct failed? "
                << flatOutput(list1.ranges()) << nl;
        }

        Info<< "unsorted: ";
        ranges.writeList(Info) << nl;

        ranges.sort();
        Info<< "sorted: ";
        ranges.writeList(Info) << nl;

        Info<< nl
            << "list linear length = " << ranges.totalSize() << nl;

        Info<< "list labels = ";
        ranges.labels().writeList(Info) << nl;

        Info<< nl;
        for (int i : { -1, 0, 5, 8, 10, 20, 26 })
        {
            Info<< "value at [" << i << "] = " << ranges[i] << nl;
        }
    }

    {
        Info<< "test intersections" << endl;
        labelRange range1(-15, 25);
        labelRange range2(7, 8);
        labelRange range3(-20, 8);
        labelRange range4(50, 8);

        Info<<range1 << " & " << range2
            << " = " << range1.subset(range2) << nl;

        Info<< range1 << " & " << range3
            << " = " << range1.subset(range3) << nl;

        Info<< range2 << " & " << range4
            << " = " << range2.subset(range4) << nl;
    }

    labelRange range;
    labelRanges ranges;

    bool removeMode = false;
    for (label argI=1; argI < args.size()-1; ++argI)
    {
        if (args[argI] == "add")
        {
            removeMode = false;
            continue;
        }
        else if (args[argI] == "del")
        {
            removeMode = true;
            continue;
        }

        {
            label start = args.get<label>(argI);
            label size  = args.get<label>(argI+1);
            ++argI;

            range.reset(start, size);
        }

        Info<< "---------------" << nl;
        if (removeMode)
        {
            Info<< "del " << range << " :";
            for (auto i : range)
            {
                Info<< " " << i;
            }
            Info<< nl;

            ranges.remove(range);
        }
        else
        {
            Info<< "add " << range  << " :";
            for (auto i : range)
            {
                Info<< " " << i;
            }
            Info<< nl;

            ranges.add(range);
        }

        Info<< "<list>" << ranges << "</list>" << nl
            << "content:";
        for (auto i : ranges)
        {
            Info<< " " << i;
        }
        Info<< nl;
    }


    {
        range.reset(5, 5);
        Info<< nl << "Tests on " << range << nl;

        printInfo(range);

        range += 3;
        Info<< "increase size " << range << nl;

        range -= 3;  // Probably not a great idea
        Info<< "decrease size " << range << nl;

        auto iter = range.begin();

        Info<< "iter: " << *iter << nl;
        ++iter;

        Info<< "iter: " << *iter << nl;

        iter += 3;
        Info<< "iter: " << *iter << nl;

        Info<< nl << "reversed:" << nl;
        forAllConstReverseIters(range, iter)
        {
            Info<< ' ' << *iter;
        }
        Info<< nl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
