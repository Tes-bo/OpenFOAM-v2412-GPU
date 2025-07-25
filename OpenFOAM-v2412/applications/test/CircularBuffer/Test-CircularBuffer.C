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
    Test-CircularBuffer

Description
    Basic tests for CircularBuffer behaviour and characteristics

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "ListOps.H"
#include "CircularBuffer.H"
#include "FlatOutput.H"

using namespace Foam;

template<class T>
inline Ostream& report
(
    const CircularBuffer<T>& buf,
    bool debugOutput = true
)
{
    buf.writeList(Info, 0);
    if (debugOutput)
    {
        Info<< " : ";
        buf.info(Info);
    }
    else
    {
        Info<< nl;
    }
    return Info;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    CircularBuffer<label> buf1(1); report(buf1);
    buf1.push_back(10); report(buf1);

    Info<< buf1.range_one() << nl;

    buf1.push_back(20); report(buf1);
    buf1.push_back(30); report(buf1);
    buf1.push_back(40); report(buf1);
    buf1.push_front(-50); report(buf1);
    buf1.push_back(60); report(buf1);
    buf1.push_back(labelList({70,80,90})); report(buf1);

    Info<< nl << "access: " << buf1 << nl;

    Info<< buf1[-12] << nl;

    Info<< "contains: " << buf1.contains(40) << nl;
    buf1.push_uniq(100); report(buf1);

    buf1 = Zero; report(buf1);

    buf1 = 500; report(buf1);

    while (buf1.size() > 2)
    {
        buf1.pop_front();
    }
    report(buf1);

    buf1.push_back(identity(5));
    buf1.emplace_front(-1000);
    buf1.emplace_back(1000);
    report(buf1);

    buf1.info(Info);
    Info<< buf1 << nl;

    CircularBuffer<label> buf2(15);
    report(buf2);

    buf2 = std::move(buf1);
    Info<< "buf1: "; report(buf1);
    Info<< "buf2: "; report(buf2);

    Info<< "for-range:";
    for (const label val : buf2)
    {
        Info<< ' ' << val;
    }
    Info<< endl;

    {
        auto iter = buf2.cbegin();
        auto endIter = buf2.cend();

        Info<< "iterated:";
        while (iter != endIter)
        {
            Info<< ' ' << *(++iter);
        }
        Info<< endl;
    }

    Info<< nl << "list: " << flatOutput(buf2.list()) << nl;

    Info<< "normal: " << flatOutput(buf2) << nl;
    buf2.reverse();
    Info<< "reverse: " << flatOutput(buf2) << nl;

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
