/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "IOstreams.H"
#include "ISLList.H"
#include "List.H"
#include "FlatOutput.H"
#include "ListOps.H"
#include "OSspecific.H"

using namespace Foam;

class Scalar
:
    public ISLList<Scalar>::link
{
public:

    // static bool verbose;

    scalar data_;

    Scalar() : data_(0) {}
    Scalar(scalar val) : data_(val) {}

    // ~Scalar() {}

    scalar value() const noexcept { return data_; }
    scalar& value() noexcept { return data_; }

    friend Ostream& operator<<(Ostream& os, const Scalar& item)
    {
        os  << item.value();
        return os;
    }
};

// bool Scalar::verbose = true;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    ISLList<Scalar> myList(new Scalar(0));

    for (int i = 0; i<10; i++)
    {
        myList.push_back(new Scalar(1.3*i));
    }

    myList.push_back(new Scalar(100.3));
    myList.push_back(new Scalar(500.3));

    Info<< "ISLList<scalar>" << myList << nl;
    Info<< nl << "flat-output: " << flatOutput(myList) << nl;

    Info<< nl << "range-for:" << nl;
    for (const auto& val : myList)
    {
        Info<< "  " << val << nl;
        // Info<<" is " << typeid(val).name() << endl;
    }

    Info<< nl << "const_iterator:" << nl;

    const ISLList<Scalar>& const_myList = myList;

    forAllConstIters(const_myList, iter)
    {
        Info<< "  " << *iter << endl;
    }

    {
        Info<< nl << "Remove element:" << nl;

        Scalar *iter = myList.removeHead();

        Info<< "  remove " << *iter;
        Info<< " => " << flatOutput(myList) << nl;

        delete iter;
    }


    Info<< nl << "Transfer: " << nl;
    Info<< "original: " << flatOutput(myList) << endl;

    ISLList<Scalar> newList;
    newList.transfer(myList);

    Info<< nl
        << "source: " << flatOutput(myList) << nl
        << "target: " << flatOutput(newList) << endl;

    myList.swap(newList);

    Info<< nl << "swap: " << nl;
    Info<< nl
        << "source: " << flatOutput(myList) << nl
        << "target: " << flatOutput(newList) << endl;

    myList.swap(newList);

    Info<< nl << "Move Construct: " << nl;

    ISLList<Scalar> list2(std::move(newList));

    Info<< nl
        << "in : " << flatOutput(newList) << nl
        << "out: " << flatOutput(list2) << nl;

    // Move back
    Info<< nl << "Move Assignment: " << nl;

    newList = std::move(list2);

    Info<< nl << "move assign: " << nl;
    Info<< nl
        << "source: " << flatOutput(list2) << nl
        << "target: " << flatOutput(newList) << endl;


    Info<< nl << "Bye." << endl;
    return 0;
}


// ************************************************************************* //
