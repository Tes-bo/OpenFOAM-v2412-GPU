/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-UList

Description
    Simple tests for UList constructors

See also
    Foam::List

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "IOstreams.H"
#include "SpanStream.H"

#include "labelList.H"
#include "ListOps.H"
#include "SubList.H"
#include "FlatOutput.H"

using namespace Foam;

template<class ListType>
void print(const ListType& list)
{
    Info << flatOutput(list) << " data addr: " << name(list.cdata()) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    List<label> source = identity(7);
    List<label> other = identity(7);

    // Text for reading as a SLList"
    string inputSLList("(10 20 30 40 50 60 70)");

    // Text for reading as a SLList"
    string inputCompound("List<label> (-1 -2 -3 -4 -5 -6 -7)");

    reverse(other);

    UList<label> ulist(source.data(), source.size());

    Info<<"source: "; print(source);
    Info<<"other:  "; print(other);
    Info<<"UList: "; print(ulist);

    {
        Info<<"shallow copy" << nl;
        ulist.shallowCopy(other);

        Info<<"source: "; print(source);
        Info<<"other:  "; print(other);
        Info<<"UList: "; print(ulist);
    }

    {
        Info<<"deep copy" << nl;
        ulist.deepCopy(source);

        Info<<"source: "; print(source);
        Info<<"other:  "; print(other);
        Info<<"UList: "; print(ulist);
    }

    {
        Info<<"Read from " << inputSLList << nl;

        ISpanStream is(inputSLList);
        is >> ulist;

//         Info<<"source: "; print(source);
//         Info<<"other:  "; print(other);
        Info<<"UList: "; print(ulist);
    }


    {
        Info<<"Read from " << inputCompound << nl;

        ISpanStream is(inputCompound);
        is >> ulist;

//         Info<<"source: "; print(source);
//         Info<<"other:  "; print(other);
        Info<<"UList: "; print(ulist);
    }

    {
        Foam::identityOp ident;

        Info<< nl << "identityOp as functor or array/map" << nl;

        for (const label val : labelRange(5))
        {
            Info<< "value:" << val
                << " () = " << ident(val)
                << " [] = " << ident[val] << nl;
        }
    }

    Info<< "\nEnd\n" << nl;
    return 0;
}


// ************************************************************************* //
