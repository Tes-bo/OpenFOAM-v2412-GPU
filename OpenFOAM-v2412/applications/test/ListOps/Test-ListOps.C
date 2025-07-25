/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2013 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Test-ListOps

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "SubList.H"
#include "ListOps.H"
#include "labelField.H"
#include "MinMax.H"
#include "face.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test Rotations:" << nl << endl;

    List<label> forwardRotate(identity(5));
    face testFace(identity(4));

    for (label i = 0; i < 8; ++i)
    {
        Info<< "Rotate forward by " << i << " : "
            << rotateList(forwardRotate, i) << endl;
    }

    for (label i = 0; i < 8; ++i)
    {
        Info<< "Rotate backward by " << i << " : "
            << rotateList(forwardRotate, -i) << endl;
    }

    Info<< nl << "Face                 : " << testFace << endl;
    Info<< "Rotate by 2          : " << rotateList(testFace, 2) << endl;
    inplaceRotateList<List, label>(testFace, -6);
    Info<< "Rotate inplace by -6 : " << testFace << nl << endl;

    Info<< "Test inplace rotate      : " << forwardRotate << endl;
    inplaceRotateList(forwardRotate, 2);
    Info<< "Rotate to the right by 2 : " << forwardRotate << endl;
    inplaceRotateList(forwardRotate, -2);
    Info<< "Rotate to the left by 2  : " << forwardRotate << endl;

    List<label> subRotate(identity(10));
    SubList<label> subL(subRotate, 5, 3);

    Info<< "Test inplace rotate on sublist : " << subRotate << endl;
    inplaceRotateList(subL, 3);
    Info<< "Rotate to the right by 3       : " << subRotate << endl;
    inplaceRotateList(subL, -8);
    Info<< "Rotate to the left by 3        : " << subRotate << endl;

    Info<< nl << nl << "Test Reversing:" << nl << endl;

    Info<< "List    : " << identity(5) << endl;
    Info<< "Reverse : " << reverseList(identity(5)) << endl;
    Info<< "List    : " << identity(6) << endl;
    Info<< "Reverse : " << reverseList(identity(6)) << nl << endl;

    List<label> test1(identity(5));
    Info<< "List            : " << test1 << endl;
    inplaceReverseList(test1);
    Info<< "Inplace Reverse : " << test1 << nl << endl;

    List<label> test2(identity(6));
    Info<< "List            : " << test2 << endl;
    inplaceReverseList(test2);
    Info<< "Inplace Reverse : " << test2 << nl << endl;

    face test3(identity(6));
    Info<< "Face            : " << test3 << endl;
    inplaceReverseList(test3);
    Info<< "Inplace Reverse : " << test3 << nl << endl;

    FixedList<label, 6> test4(identity(6));
    Info<< "FixedList       : " << test4 << endl;
    inplaceReverseList(test4);
    Info<< "Inplace Reverse : " << test4 << nl << endl;

    List<label> test5(identity(9));
    SubList<label> test5SubList(test5, 4, 3);
    Info<< "List                            : " << test5 << endl;
    inplaceReverseList(test5SubList);
    Info<< "Reverse Sublist between 3 and 6 : " << test5 << nl << endl;

    Info<< nl << "Test lambda predicates:" << nl << endl;

    List<label> test6(identity(11, -4));

    // Add multiplier for mpore interesting testing
    std::for_each(test6.begin(), test6.end(), [](label& x){ x *= 3; });

    // Randomize the list
    Foam::shuffle(test6);

    Info<< "randomized input list: " << flatOutput(test6) << nl;

    const auto evenNonZero = [](const label& x){ return x && !(x % 2); };

    Info<< "location of first even/non-zero: "
        << ListOps::find_if(test6, evenNonZero) << nl;

    Info<< "find > 12 && divisible by 5 : "
        << ListOps::find_if
           (
               test6,
               [](const label& x) { return x > 12 && !(x % 5); }
           ) << nl;

    Info<< "Found >= 8 : "
        << ListOps::found_if(test6, labelMinMax(8, labelMax)) << nl;

    Info<< "Found >= 25 : "
        << ListOps::found_if(test6, labelMinMax(25, labelMax)) << nl;


    Info<< "Subset of non-zero, even values: "
        << subsetList(test6, evenNonZero) << nl
        << endl;

    Info<< "Subset of non-zero, even values: "
        << IndirectList<label>::subset_if(test6, evenNonZero) << nl
        << endl;

    test6.push_back(identity(13, 12));

    Info<< "Randomized: " << flatOutput(test6) << endl;
    inplaceUniqueSort(test6);
    Info<< "Unique    : " << flatOutput(test6) << endl;


    // List reorder
    labelList oldToNew(identity(40));
    Foam::shuffle(oldToNew);

    // Force a few -1:
    oldToNew[4] = oldToNew[8] = -1;

    Info<<"Test reorder - oldToNew:" << nl
        << flatOutput(oldToNew) << nl << nl;

    bitSet bitset
    (
        ListOps::createWithValue<bool>
        (
            25,
            labelList({8, 12, 15, 22, 4}),
            true
        )
    );

    Info<<"bitset input: " << flatOutput(bitset) << nl;
    inplaceReorder(oldToNew, bitset);
    Info<<"     reorder: " << flatOutput(bitset) << nl << nl;

    PackedList<2> packed
    (
        ListOps::createWithValue<label>
        (
            25,
            labelList({8, 12, 15, 22, 4}),
            2
        )
    );

    Info<< "packed input: " << packed << nl;
    inplaceReorder(oldToNew, packed);
    Info<<"     reorder: " << packed << nl << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
