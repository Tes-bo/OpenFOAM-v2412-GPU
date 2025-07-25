/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Test-Hashing2

Description

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

#include "Hash.H"
#include "HashSet.H"

using namespace Foam;

void infoHashString
(
    unsigned modulus,
    std::initializer_list<std::string> lst
)
{
    if (modulus)
    {
        Info<< "basic string hashing (mod " << label(modulus) << ")" << nl;

        for (const auto& str : lst)
        {
            Info<<"hash(" << str.c_str() << ")="
                << (Hash<string>()(str) % modulus) << nl;
        }

    }
    else
    {
        Info<< "basic string hashing" << nl;

        for (const auto& str : lst)
        {
            Info<<"hash(" << str.c_str() << ")="
                << Hash<string>()(str) << nl;
        }
    }
}



void reportHashList(const UList<string>& list)
{
    Info<< "contiguous = " << is_contiguous<string>::value << nl << nl;

    for (const string& val : list)
    {
        unsigned hash1 = string::hasher()(val);

        Info<< hex << hash1 << ": " << val << nl;
    }
}


void reportHashList(const UList<label>& list)
{
    Info<<"contiguous = " << is_contiguous<label>::value << nl << nl;

    for (const label val : list)
    {
        // Direct value
        unsigned hash1 = Hash<label>()(val);

        // Hashed byte-wise
        unsigned hash2 = Hash<label>()(val, 0);

        Info<< hex << hash1
            << " (seeded: " << hash2 << ")"
            << ": " << dec << val << nl;
    }
}


void reportHashList(const UList<face>& list)
{
    Info<<"contiguous = " << is_contiguous<label>::value << nl << nl;

    for (const face& f : list)
    {
        // Direct value
        unsigned hash1 = face::hasher()(f);

        unsigned hash2 = Hash<face>()(f);

        Info<< hex << "face hash " << hash1
            << " Hash<face> " << hash2
            << ": " << dec << flatOutput(f) << nl;
    }
}


void reportHashList(const UList<labelList>& list)
{
    for (const labelList& val : list)
    {
        unsigned hash1 = Foam::Hasher(val.cdata(), val.size_bytes());

        unsigned hash2 = Hash<labelList>()(val);
        unsigned hash2b = labelList::hasher()(val);

        Info<< hex << hash1 << " or " << hash2
            << "(" << hash2b << ") "
            << ": " << dec << val << nl;
    }

    unsigned hash2 = Hash<labelListList>()(list);
    unsigned hash2bad = Foam::Hasher(&list, sizeof(list));

    Info<< hex << hash2 << " : " << dec << flatOutput(list) << nl
        << hex << hash2bad << " as direct hash would be wrong"
        << dec << nl;
}


void reportHashList(const UList<wordPair>& list)
{
    Info<<"contiguous = " << is_contiguous<wordPair>::value << nl << nl;

    for (const wordPair& pr : list)
    {
        unsigned hash1 = Hash<wordPair>()(pr);

        // as FixedList
        unsigned hash2 = wordPair::hasher()(pr);

        // as FixedList
        unsigned hash2sym = wordPair::symmHasher()(pr);

        // as FixedList
        unsigned hash3 = Hash<FixedList<word,2>>()(pr);

        Info<< hex << hash1 << " (as FixedList: " << hash2
            << ") or " << hash3
            << " symm-hash:" << hash2sym
            << " : "<< dec << flatOutput(pr) << nl;
    }
}


void reportHashList(const UList<labelPair>& list)
{
    Info<<"contiguous = " << is_contiguous<labelPair>::value << nl << nl;

    for (const labelPair& pr : list)
    {
        unsigned hash1 = Hash<labelPair>()(pr);

        // as FixedList
        unsigned hash2 = labelPair::hasher()(pr);

        // as FixedList
        unsigned hash3 = Hash<labelPair>()(pr);

        Info<< hex << hash1 << " (as FixedList: " << hash2
            << ") or " << hash3
            << " : "<< dec << pr << nl;
    }
}


void reportHashList(const UList<labelPairPair>& list)
{
    Info<<"contiguous = " << is_contiguous<labelPairPair>::value << nl << nl;

    for (const labelPairPair& pr : list)
    {
        unsigned hash1 = Hash<labelPairPair>()(pr);

        // as FixedList
        unsigned hash2 = labelPairPair::hasher()(pr);

        // as FixedList
        unsigned hash3 = Hash<labelPairPair>()(pr);

        Info<< hex << hash1 << " (as FixedList: " << hash2
            << ") or " << hash3
            << " : "<< dec << pr << nl;
    }
}


void reportHashList(const UList<edge>& list)
{
    Info<<"contiguous = " << is_contiguous<edge>::value << nl << nl;

    for (const edge& e : list)
    {
        unsigned hash1 = Hash<edge>()(e);

        // as FixedList
        unsigned hash2 = labelPair::hasher()(e);

        // as FixedList
        unsigned hash3 = Hash<labelPair>()(e);

        Info<< hex << hash1 << " (as FixedList: " << hash2
            << ") or " << hash3
            << " : "<< dec << e << nl;
    }
}


void reportHashList(const UList<triFace>& list)
{
    Info<<"contiguous = " << is_contiguous<triFace>::value << nl << nl;

    for (const triFace& f : list)
    {
        // direct value
        unsigned hash1 = Hash<triFace>()(f);
        unsigned hash2 = FixedList<label, 3>::hasher()(f);

        Info<< hex << hash1 << " (as FixedList: " << hash2
            << "): " << dec << f << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    infoHashString(8, {"asdathis1", "adsxf", "hij", "klmpq"});

    IFstream is("hashingTests");

    if (is.good())
    {
        Info<< nl << "Process " << is.name() << " file ..." << nl;
    }
    else
    {
        Info<< nl << "No " << is.name() << " file found ..." << nl;
    }

    token tok;

    while (is.good() && is.read(tok) && tok.good())
    {
        const word listType(tok.wordToken());

        Info<< nl;
        IOobject::writeDivider(Info) << listType << nl;

        if (listType == "stringList")
        {
            stringList list(is);
            reportHashList(list);
        }
        else if (listType == "labelList")
        {
            labelList list(is);
            reportHashList(list);
        }
        else if (listType == "faceList")
        {
            faceList list(is);
            reportHashList(list);
        }
        else if (listType == "labelListList")
        {
            List<labelList> list(is);
            reportHashList(list);
        }
        else if (listType == "wordPairList")
        {
            List<wordPair> list(is);
            reportHashList(list);
        }
        else if (listType == "labelPairList")
        {
            labelPairList list(is);
            reportHashList(list);
        }
        else if (listType == "labelPairPairList")
        {
            List<labelPairPair> list(is);
            reportHashList(list);
        }
        else if (listType == "edgeList")
        {
            edgeList list(is);
            reportHashList(list);
        }
        else if (listType == "triFaceList")
        {
            triFaceList list(is);
            reportHashList(list);
        }
        else
        {
            Info<< "unknown type: " << listType << nl;
        }
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
