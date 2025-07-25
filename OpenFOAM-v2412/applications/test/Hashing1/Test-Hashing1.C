/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test-Hashing1

Description
    Test/verify overloads of Hash function

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"

#include "stringList.H"
#include "labelList.H"
#include "labelPair.H"
#include "wordPair.H"
#include "edgeList.H"
#include "faceList.H"
#include "triFaceList.H"

#define FULLDEBUG
#include "HashFunction.H"

using namespace Foam;

template<class T>
unsigned rawHasher(const T& obj)
{
    return Foam::Hasher(&obj, sizeof(T), 0u);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        typedef unsigned Type;
        Type value = 100;

        Info<< "hash " << typeid(value).name() << " of " << value << nl;
        Info<< "    Hasher: " << rawHasher(value) << nl;
        Info<< "    Hash<>: " << Hash<Type>()(value) << nl;
    }
    {
        typedef int32_t Type;
        Type value = 100;

        Info<< "hash " << typeid(value).name() << " of " << value << nl;
        Info<< "    Hasher: " << rawHasher(value) << nl;
        Info<< "    Hash<>: " << Hash<Type>()(value) << nl;
    }
    {
        typedef int64_t Type;
        Type value = 100;

        Info<< "hash " << typeid(value).name() << " of " << value << nl;
        Info<< "    Hasher: " << rawHasher(value) << nl;
        Info<< "    Hash<>: " << Hash<Type>()(value) << nl;
    }

    HashFun<std::string>().info();
    HashFun<string>().info();
    HashFun<Foam::word>().info();
    HashFun<Foam::keyType>().info();

    HashFun<int>().info();
    HashFun<label>().info();
    HashFun<float>().info();
    HashFun<double>().info();

    {
        float value = 15.f;

        Info<< "hash of " << Foam::name(&value)
            << " = " << HashFun<void*>()(&value) << nl;
    }

    HashFun<labelList>().info();
    HashFun<wordUList>().info();
    HashFun<edge>().info();

    HashFun<Pair<label>>().info();
    HashFun<labelPair>().info();
    HashFun<labelPairPair>().info();

    HashFun<Tuple2<label, word>>().info();

    {
        typedef Tuple2<label, word> Type;
        Type obj(10, "test");
        Info<< obj << " hash=" << Hash<Type>()(obj) << nl;
    }

    {
        typedef Tuple2<label, label> Type;
        Type obj(10, 12);
        Info<< obj << " hash=" << Hash<Type>()(obj) << nl;
    }

    {
        typedef Pair<label> Type;
        Type obj(10, 12);
        Info<< obj << " hash=" << Hash<Type>()(obj) << nl;
    }

    {
        typedef std::pair<label, label> Type;
        Type obj(10, 12);
        Info<< obj << " hash=" << Hash<Type>()(obj) << nl;

        HashSet<Type> hs;
        hs.insert(obj);
        hs.erase(obj);
    }

    {
        Pair<label>::hasher op;
        Info<< "hasher: " << op(Pair<label>(10, 12)) << nl;
    }

    // Not supported
    #if 0
    {
        Tuple2<label, label>::hasher op;
        Info<< "hasher: " << op(Tuple2<label>(10, 12)) << nl;
    }
    #endif


    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
