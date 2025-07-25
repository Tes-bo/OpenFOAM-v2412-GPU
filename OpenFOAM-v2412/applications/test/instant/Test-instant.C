/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

Description
    Test instant, fileNameInstant

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "instant.H"
#include "Pair.H"
#include "fileNameInstant.H"
#include "DynamicList.H"

using namespace Foam;

template<class T>
Ostream& printInstant(const UList<T>& times, const label i)
{
    if (i >= 0 && i < times.size())
    {
        Info<< " (" << times[i] << ")";
    }
    return Info;
}

template<class T>
Ostream& printInstant(const UList<T>& times, const Pair<label>& range)
{
    printInstant(times, range.first());
    printInstant(times, range.second());
    return Info;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    DynamicList<instant> times;

    times.emplace_back();
    times.emplace_back(12, "abc");
    times.emplace_back(3.14159);
    times.emplace_back(300.456, "def");
    times.emplace_back(454.456, "xyz");
    times.emplace_back(10, "ten");
    times.push_back({15, "fifteen"});

    {
        word timeName("twenty");

        Info<<"move append: <" << timeName << '>' << nl;
        times.emplace_back(20, std::move(timeName));
        Info<<"after append: <" << timeName << '>' << nl;
    }

    Info<< nl << "times:" << times << nl;
    sort(times);
    Info<< "Sorted:" << times << nl;

    for (const scalar val : { -0.5, 5.0, 18.0, 25.0, 450.0, 480.0 })
    {
        label start = instant::findStart(times, val);
        Pair<label> range = instant::findRange(times, val);

        Info<< nl
            << "time:" << val << nl;
        Info<< "    start:" << start;
        printInstant(times, start) << nl;

        Info<< "    range:" << range;
        printInstant(times, range) << nl;
    }

    DynamicList<fileNameInstant> files;
    files.emplace_back();
    files.emplace_back(12, "twelve");
    files.emplace_back(3.14, "/path/almost-pi");
    files.emplace_back(300, "/dev/value");
    files.emplace_back(454, "/tmp/xyz");
    files.emplace_back(10, "ten");

    Info<< nl << "files:" << files << nl;
    sort(files);
    Info<< "Sorted:" << files << nl;


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
