/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-stdFoam-span

Description
    Basic functionality test for span

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "IOstreams.H"
#include "scalar.H"
#include "vector.H"
#include "labelRange.H"
#include "scalarList.H"

using namespace Foam;

template<class Type>
void printInfo(const UList<Type>& list)
{
    Info<< "list: " << flatOutput(list) << nl
        << "data: " << Foam::name(list.cdata())
        << " size: " << Foam::name(list.size()) << nl;
}


template<class Type>
void printInfo(const stdFoam::span<Type> span)
{
    Info<< "span: " << Foam::name(span.data())
        << ", " << span.size() << nl;
}

template<class Type>
void printContent(const stdFoam::span<Type> span)
{
    Info<< "range:";
    for (const auto& val : span)
    {
        Info<< " " << val;
    }
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"

    {
        const List<label> list1(identity(15));

        printInfo(list1);

        stdFoam::span<const label> span1;
        printInfo(span1);

        span1 = stdFoam::span<const label>(list1.cdata(), list1.size());
        printInfo(span1);
        printContent(span1);

        auto subspan1 = span1.subspan(5);

        printInfo(subspan1);
        printContent(subspan1);

        subspan1 = span1.subspan(5, 30);
        printInfo(subspan1);
        printContent(subspan1);

        // non-const access: span1[6] = -100;
        // non-const access: Info<< Foam::name(span1.data_bytes()) << endl;

        Info<< Foam::name(span1.cdata_bytes()) << endl;
    }

    {
        List<label> list1(identity(15));

        printInfo(list1);

        stdFoam::span<label> span1;
        printInfo(span1);

        span1 = stdFoam::span<label>(list1.data(), list1.size());
        printInfo(span1);
        printContent(span1);

        auto subspan1 = span1.subspan(5);

        printInfo(subspan1);
        printContent(subspan1);

        subspan1 = span1.subspan(5, 30);
        printInfo(subspan1);
        printContent(subspan1);

        span1[6] = -100;
        printInfo(span1);
        printContent(span1);

        Info<< Foam::name(span1.data_bytes()) << endl;


        // With const span
        const auto span2 = stdFoam::span<label>(list1.data(), list1.size());
        printContent(span2);

        // Span may be const, but its contents (pointers) needn't be
        span2[1] = 100;
        span2.front() = -1000;
        span2.back() = 1000;
        printInfo(list1);

        printContent(span2.first(5));
        printContent(span2.last(5));
    }

    Info<< "\nDone\n" << nl;
    return 0;
}

// ************************************************************************* //
