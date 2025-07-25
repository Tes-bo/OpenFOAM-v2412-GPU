/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-parallel-broadcast

Description
    Test for various broadcast routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "bitSet.H"
#include "vector.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

using namespace Foam;


template<class T>
void printPre(const T& value)
{
    Info<< nl << "is_contiguous:" << is_contiguous<T>::value << endl;
    Pout<< "pre-broadcast: " << value << endl;
}

template<class T>
void printPost(const T& value)
{
    Pout<< "post-broadcast: " << value << endl;
}


template<class T>
void testBroadcast(T& value)
{
    printPre(value);
    Pstream::broadcast(value);
    printPost(value);
}

template<class T>
void testBroadcast(List<T>& values)
{
    Info<< nl << "is_contiguous:" << is_contiguous<T>::value << endl;
    Pout<< "pre-broadcast: " << flatOutput(values) << endl;
    Pstream::broadcast(values);
    Pout<< "post-broadcast: " << flatOutput(values) << endl;
}


void testBroadcast(bitSet& values)
{
    Pout<< "pre-broadcast: "
        << values.size() << ": " << flatOutput(values.values()) << endl;
    Pstream::broadcast(values);
    Pout<< "post-broadcast: "
        << values.size() << ": " << flatOutput(values.values()) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    #include "setRootCase.H"
    #include "createTime.H"

    {
        label value = -1;
        if (Pstream::master())
        {
            value = UPstream::nProcs();
        }
        testBroadcast(value);
    }

    {
        labelList values;
        if (Pstream::master())
        {
            values = identity(UPstream::nProcs());
        }
        testBroadcast(values);
    }

    {
        word value;
        if (Pstream::master())
        {
            value = args.executable();
        }
        printPre(value);
        Pstream::broadcast(value);   // Streamed broadcast
        printPost(value);
    }

    {
        wordList values;
        if (Pstream::master())
        {
            values.resize(UPstream::nProcs());
            forAll(values, i)
            {
                values[i] = "value_" + Foam::name(i);
            }
        }
        testBroadcast(values);
    }

    {
        vector values(vector::uniform(-1));
        if (Pstream::master())
        {
            values = vector(1,2,3);
        }
        testBroadcast(values);
    }

    {
        FixedList<vector, 3> values(vector::uniform(-1));
        if (Pstream::master())
        {
            values = vector(1,2,3);

            scalar mult = 1;
            for (auto& v : values)
            {
                v *= mult;
                mult += 1;
            }
        }
        testBroadcast(values);
    }

    {
        bitSet values;
        if (Pstream::master())
        {
            values.set(labelList({1, 4, 8}));
            values.resize(10);
        }
        else
        {
            // Just something different
            values.set(labelList({0, 2}));
            values.resize(5);
        }
        testBroadcast(values);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
