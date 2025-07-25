/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Test-vectorSpeed

Description
    Test speeds, usability of some field operations

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "primitiveFields.H"
#include "cpuTime.H"
#include "IOstreams.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("lerp");

    argList args(argc, argv);

    const label nIter = 1000;
    const label size = (1000000);

    Info<< "Initialising fields. size:" << size
        << " max:" << labelMax << endl;

    scalarField onet(size);

    vectorField
        vf1(size, vector::one),
        vf2(size, vector::one),
        vf3(size, vector::one),
        vf4(size);

    Info<< "Start loop: " << nIter << endl;

    cpuTime timing;

    // Timing is mostly malloc anyhow...

    if (!args.found("lerp"))
    {
        Info<< "vectorField algebra" << endl;

        for (int j=0; j<nIter; j++)
        {
            vf4 = vf1 + vf2 - vf3;
        }

        Info<< "Timing = " << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;
    }

    if (args.found("lerp"))
    {
        Info<< nl
            << "testing timings with lerp 0.75" << nl
            << "Main bottlenecks: malloc and looping" << nl;

        scalar t = 0.75;
        scalarField fract(size, t);

        // Basic compilation tests:
        lerp(vf1, vf2, fract)().size();             // list list list
        lerp(0.5*vf1, vf2, fract)().size();         // tmp  list list
        lerp(vf1, 0.5*vf2, fract)().size();         // list tmp  list
        lerp(vf1, vf2, 0.5*fract)().size();         // list list tmp
        lerp(0.5*vf1, 0.5*vf2, fract)().size();     // tmp  tmp  list
        lerp(0.5*vf1, vf2, 0.5*fract)().size();     // tmp  list tmp
        lerp(vf1, 1.0*vf2, 0.5*fract)().size();     // list tmp  tmp
        lerp(0.5*vf1, 0.5*vf2, 0.5*fract)().size(); // tmp  tmp  tmp

        // Basic compilation tests:
        lerp(vf1, vf2, t)().size();                 // list list scalar
        lerp(0.5*vf1, vf2, t)().size();             // tmp  list scalar
        lerp(vf1, 0.5*vf2, t)().size();             // list tmp  scalar
        lerp(0.5*vf1, 0.5*vf2, t)().size();         // tmp  tmp  scalar
        lerp(vf1, 1.0*vf2, 0.5*t)().size();         // list tmp  scalar

        // Other combinations (not yet included)
        // list scalar list...
        // scalar list list...


        // Probable base line of single instruction multiple data
        for (int j=0; j<nIter; j++)
        {
            const label loopLen = (vf1).size();

            for (label i = 0; i < loopLen; ++i)
            {
                vf4[i] = (1-fract[i])*vf1[i] + fract[i]*vf2[i];
            }
        }

        Info<< "(1-t)*a + t*b : Timing (no malloc) = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;

        // With all loops written out
        // - corresponds to what we would normally have
        for (int j=0; j<nIter; j++)
        {
            const label loopLen = (vf1).size();

            for (label i = 0; i < loopLen; ++i)
            {
                onet[i] = (1-fract[i]);
            }

            for (label i = 0; i < loopLen; ++i)
            {
                vf3[i] = (fract[i] * vf2[i]);
            }

            for (label i = 0; i < loopLen; ++i)
            {
                vf4[i] = onet[i] * vf1[i];
            }

            for (label i = 0; i < loopLen; ++i)
            {
                vf4[i] += vf3[i];
            }
        }

        Info<< "(1-t)*a + t*b : Looping (no malloc) = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;

        for (int j=0; j<nIter; j++)
        {
            lerp(vf4, vf1, vf2, fract);
        }

        Info<< "lerp(a, b, t) : Timing (no malloc) = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;

        Info<< nl;

        // Now with tmp fields
        // ~~~~~~~~~~~~~~~~~~~

        for (int j=0; j<nIter; j++)
        {
            vf4 = (1-fract)*vf1 + fract*vf2;
        }

        Info<< "(1-t)*a + t*b : Timing (with malloc) = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;

        for (int j=0; j<nIter; j++)
        {
            // With field 't'
            vf4 = lerp(vf1, vf2, fract);
        }

        Info<< "lerp(a, b, t) : Timing (with malloc) = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;

        for (int j=0; j<nIter; j++)
        {
            // With scalar 't'
            vf4 = lerp(vf1, vf2, t);
        }

        Info<< "lerp(a, b, t) : Timing (with malloc) scalar = "
            << timing.cpuTimeIncrement() << " s" << endl;

        Snull<< vf4[1] << endl << endl;
    }


    Info<< nl << "Done" << nl << endl;
    return 0;
}


// ************************************************************************* //
