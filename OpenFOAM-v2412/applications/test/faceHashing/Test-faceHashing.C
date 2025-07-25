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
    Test-faceHashing

Description
    Basic tests of face/triFace hashing

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"

#include "faceList.H"
#include "triFaceList.H"
#include "HashSet.H"
#include "HashTable.H"

using namespace Foam;

template<class HashSetType, class FaceListType>
void checkHashSet(const HashSetType& hs, const FaceListType& faces)
{
    Info<< hs << nl;
    for (const auto& f : faces)
    {
        Info<< "  " << f << " found=" << hs.found(f) << nl;
    }
    Info<< nl;
}


void checkHashes(const faceList& faces)
{
    face::hasher op1;
    face::symmHasher op2;

    Info<< "face hasher/symmHasher: " << faces.size() << " faces" << nl;
    for (const face& f : faces)
    {
        Info<< "  " << f << " symmhash=" << op2(f)
            << " hash=" << op1(f) << nl;
    }
    Info<< nl;
}


void checkHashes(const triFaceList& faces)
{
    triFace::hasher op1;

    Info<< "triFace hasher: " << faces.size() << " faces" << nl;
    for (const triFace& f : faces)
    {
        Info<< "  " << f << " hash=" << op1(f) << nl;
    }
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        faceList faces
        ({
            face{0, 1, 2, 3},  // regular
            face{0, 3, 2, 1},  // flip
            face{2, 3, 0, 1},  // rotate
            face{2, 1, 0, 3},  // rotate (flip)
            face{2, 0, 1, 3},  // permute
            face{2, 1, 0, 3},  // permute
        });

        checkHashes(faces);

        HashSet<face, face::hasher> hs1;
        HashSet<face, face::symmHasher> hs2;
        hs1.insert(faces);
        hs2.insert(faces);

        Info<< "hashset (hasher)" << nl;
        checkHashSet(hs1, faces);
        Info<< nl;

        hs1.erase(faces[0]);
        Info<< "remove " << faces[0] << nl;
        Info<< "hashset (hasher)" << nl;
        checkHashSet(hs1, faces);
        Info<< nl;

        Info<< "hashset (symmHasher)" << nl;
        checkHashSet(hs2, faces);
        Info<< nl;

        hs2.erase(faces[0]);
        Info<< "remove " << faces[0] << nl;
        Info<< "hashset (symmHasher)" << nl;
        checkHashSet(hs2, faces);
        Info<< nl;
    }

    {
        triFaceList faces
        ({
            triFace{0, 1, 2},  // regular
            triFace{0, 2, 1},  // flip
            triFace{2, 0, 1},  // rotate
        });

        checkHashes(faces);

        HashSet<triFace> hs1;
        hs1.insert(faces);

        Info<< "hashset (symmHasher)" << nl;
        checkHashSet(hs1, faces);
        Info<< nl;
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
