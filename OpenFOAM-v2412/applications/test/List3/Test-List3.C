/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Test-List3

Description
    Test list construction

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

template<class T>
void printAddress(const UList<T>& list)
{
    Info<< "list addr: " << name(&list)
        << " data addr: " << name(list.cdata()) << nl;
}


template<class T>
void printAddress(const SLList<T>& list)
{
    Info<< "list addr: " << name(&list)
        << " data addr: ???" << nl;
}


template<class T>
void printAddresses(const List<List<T>>& list)
{
    for (const auto& elem : list)
    {
        printAddress(elem);
    }
}


template<class T>
void printAddresses(const SLList<List<T>>& list)
{
    for (const auto& elem : list)
    {
        printAddress(elem);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("labelListList");

    argList args(argc, argv, false);

    if (args.options().empty())
    {
        Info<< nl << "Specify an option! " << nl << endl;
    }

    if (args.found("labelListList"))
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            if (true)
            {
                IFstream is(args.get<fileName>(argi));

                Info<< nl << nl
                    << "read from " << is.name() << nl << endl;

                SLList<List<label>> sll(is);
                Info<< "read " << sll.size() << " entries" << nl;

                Info<< "sll" << nl;
                for (const auto& elem : sll)
                {
                    printAddress(elem);
                }


//            List<List<label>> list(std::move(sll));
                List<List<label>> list;
                Info<< "move to List" << nl;
                list = std::move(sll);

                Info<< "sll" << nl;
                for (const auto& elem : sll)
                {
                    printAddress(elem);
                }
                Info<< "list" << nl;
                printAddresses(list);
            }

            if (true)
            {
                IFstream is(args.get<fileName>(argi));

                Info<< nl << nl
                    << "read from " << is.name() << nl << endl;

                List<List<label>> list(is);

                Info<< "list" << nl;
                for (const auto& elem : list)
                {
                    printAddress(elem);
                }
            }
        }
    }


    Info<< nl << "Done" << nl << endl;
    return 0;
}


// ************************************************************************* //
