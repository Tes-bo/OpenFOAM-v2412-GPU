/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Test-coordinateSystem

Description
    Expand coordinate system definitions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "coordinateSystems.H"
#include "identityRotation.H"
#include "cartesianCS.H"
#include "indirectCS.H"
#include "Fstream.H"
#include "IOstreams.H"
#include "transform.H"

using namespace Foam;


template<class T>
void testTransform(const coordinateSystem& cs, const point& p, const T& val)
{
    Info<< "    " << pTraits<T>::typeName << ": " << val
        << " transform: " << cs.transform(p, val)
        << " invTransform: " << cs.invTransform(p, val) << nl;

    // Info<< " both: " << cs.invTransform(p, cs.transform(p, val)) << nl;
}


void basicTests(const coordinateSystem& cs)
{
    cs.writeEntry(cs.name(), Info);

    if (const auto* cartptr = isA<coordSystem::cartesian>(cs))
    {
        if (!cartptr->good())
        {
            Info<< "invalid cartesian = " << (*cartptr)
                << " with: " << (*cartptr).R() << nl;
        }
    }

    Info<< "rotation: " << cs.R() << nl;

    List<point> testPoints
    ({
        {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1},
    });


    for (const point& p : testPoints)
    {
        Info<< nl
            << "  test point: " << p
            << " = local point " << cs.transformPoint(p)
            << " = local coord " << cs.localPosition(p) << nl;

        const vector v1(1, 1, 1);
        const tensor t1(tensor::I);
        const tensor t2(1, 2, 3, 4, 5, 6, 7, 8, 9);

        testTransform(cs, p, v1);
        testTransform(cs, p, t1);
        testTransform(cs, p, t2);
    }

    Info<< nl;
}


void doTest(const dictionary& dict)
{
    Info<< dict.dictName() << dict << nl;

    // Could fail?
    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOErr = FatalIOError.throwing(true);

    try
    {
        auto cs1ptr = coordinateSystem::New(dict, word::null);
        coordinateSystem& cs1 = *cs1ptr;
        cs1.rename(dict.dictName());

        basicTests(cs1);
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "Caught FatalIOError " << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOErr);
}


void doTest(const objectRegistry& obr, const dictionary& dict)
{
    Info<< dict.dictName() << dict << nl;

    // Could fail?
    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOErr = FatalIOError.throwing(true);

    try
    {
        auto cs1ptr = coordinateSystem::New(obr, dict, word::null);
        coordinateSystem& cs1 = *cs1ptr;

        basicTests(cs1);
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "Caught FatalIOError " << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOErr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    if (args.found("case"))
    {
        Info<<"using case for tests" << nl;

        #include "createTime.H"

        const coordinateSystems& systems = coordinateSystems::New(runTime);

        Info<< systems.size() << " global systems" << nl;

        for (const coordinateSystem& cs : systems)
        {
            basicTests(cs);
        }

        // systems.write();

        for (label argi=1; argi < args.size(); ++argi)
        {
            const auto dictFile = args.get<fileName>(argi);
            IFstream is(dictFile);

            dictionary inputDict(is);

            for (const entry& dEntry : inputDict)
            {
                if (dEntry.isDict())
                {
                    doTest(runTime, dEntry.dict());
                }
            }
        }
    }
    else if (args.size() <= 1)
    {
        Info<<"no coordinateSystem dictionaries to expand" << nl;
    }
    else
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            const auto dictFile = args.get<fileName>(argi);
            IFstream is(dictFile);

            dictionary inputDict(is);

            for (const entry& dEntry : inputDict)
            {
                if (dEntry.isDict())
                {
                    doTest(dEntry.dict());
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
