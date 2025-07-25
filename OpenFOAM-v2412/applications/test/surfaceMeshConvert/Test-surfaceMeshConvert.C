/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Test-surfaceMeshConvert

Group
    grpSurfaceUtilities

Description
    Test conversions from one surface mesh format to another.

Usage
    \b Test-surfaceMeshConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface

      - \par -orient
        Check face orientation on the input surface

      - \par -testModify
        Test modification mechanism

      - \par -scale \<scale\>
        Specify a scaling factor for writing the files

      - \par -triSurface
        Use triSurface library for input/output
Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "surfMesh.H"
#include "surfFields.H"

#include "MeshedSurfaces.H"
#include "ModifiableMeshedSurface.H"
#include "UnsortedMeshedSurfaces.H"
#include "SpanStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convert between surface formats, "
        "but primarily for testing functionality\n"
        "Normally use surfaceMeshConvert instead."
    );

    argList::noParallel();
    argList::addArgument("inputFile");
    argList::addArgument("outputFile");

    argList::addBoolOption
    (
        "clean",
        "perform some surface checking/cleanup on the input surface"
    );
    argList::addBoolOption
    (
        "orient",
        "check surface orientation"
    );

    argList::addBoolOption
    (
        "testModify",
        "Test modification mechanism (MeshedSurface)"
    );

    argList::addBoolOption
    (
        "surfMesh",
        "test surfMesh output"
    );
    argList::addBoolOption
    (
        "triSurface",
        "use triSurface for read/write"
    );
    argList::addBoolOption
    (
        "unsorted",
        "use UnsortedMeshedSurface instead of MeshedSurface, "
        "or unsorted output (with -triSurface option)"
    );
    argList::addBoolOption
    (
        "triFace",
        "use triFace instead of face"
    );
    argList::addBoolOption
    (
        "stdout",
        "ignore output filename and write to stdout"
    );

    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );

    #include "setRootCase.H"

    const bool     optStdout = args.found("stdout");
    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 0);

    const fileName importName = args[1];
    const fileName exportName = optStdout ? "-stdout" : args[2];

    if (importName == exportName)
    {
        FatalError
            << "Output file would overwrite input file."
            << exit(FatalError);
    }

    if
    (
        !args.found("triSurface")
    &&
         (
            !MeshedSurface<face>::canRead(importName, true)
         ||
            (
                !optStdout
             && !MeshedSurface<face>::canWriteType(exportName.ext(), true)
            )
         )
    )
    {
        return 1;
    }

    if (args.found("triSurface"))
    {
        triSurface surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area         : " << sum(surf.magSf()) << nl
            << endl;

        // check: output to ostream, construct from istream
        {
            OCharStream os;
            os << surf;
            ISpanStream is(os.view());

            // both work:
            triSurface surf2(is);

            // OR
            // is.rewind();
            // triSurface surf2;
            // is >> surf2;

            // surf2.read(is); // FAIL: private method
        }

        if (args.found("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.found("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< "Area         : " << sum(surf.magSf()) << nl
                << endl;
        }

        if (optStdout)
        {
            Info<< surf;
        }
        else
        {
            // normally write sorted (looks nicer)
            surf.write(exportName, !args.found("unsorted"));
        }
    }
    else if (args.found("unsorted"))
    {
        UnsortedMeshedSurface<face> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area         : " << sum(surf.magSf()) << nl
            << endl;

        // check: output to ostream, construct from istream
        {
            OCharStream os;
            os << surf;
            ISpanStream is(os.view());

            // both work:
            UnsortedMeshedSurface<face> surf2(is);

            // OR
            // is.rewind();
            // UnsortedMeshedSurface<face> surf2;
            // is >> surf2;

            // surf2.readIstream(is);  // FAIL: private method
        }

        if (args.found("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.found("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< "Area         : " << sum(surf.magSf()) << nl
                << endl;
        }

        if (optStdout)
        {
            Info<< surf;
        }
        else
        {
            surf.write(exportName);
        }
    }
    else if (args.found("triFace"))
    {
        MeshedSurface<triFace> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area         : " << sum(surf.magSf()) << nl
            << endl;

        // check: output to ostream, construct from istream
        {
            OCharStream os;
            os << surf;
            ISpanStream is(os.view());

            // both work:
            MeshedSurface<face> surf2(is);

            // OR
            // is.rewind();
            // MeshedSurface<face> surf2;
            // is >> surf2;

            // surf2.read(is);  // FAIL: private method
        }

        if (args.found("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.found("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< "Area         : " << sum(surf.magSf()) << nl
                << endl;
        }

        if (optStdout)
        {
            Info<< surf;
        }
        else
        {
            surf.write(exportName);
        }
    }
    else
    {
        MeshedSurface<face> surf(importName);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area         : " << sum(surf.magSf()) << nl
            << endl;

        // check: output to ostream, construct from istream
        {
            OCharStream os;
            os << surf;
            ISpanStream is(os.view());

            // both work:
            MeshedSurface<face> surf2(is);

            // OR
            // is.rewind();
            // MeshedSurface<face> surf2;
            // is >> surf2;

            // surf2.read(is);  // FAIL: private method
        }

        if (args.found("orient"))
        {
            Info<< "Checking surface orientation" << endl;
            PatchTools::checkOrientation(surf, true);
            Info<< endl;
        }

        if (args.found("clean"))
        {
            Info<< "Cleaning up surface" << endl;
            surf.cleanup(true);
            surf.writeStats(Info);
            Info<< endl;
        }

        if (args.found("testModify"))
        {
            Info<< "Use ModifiableMeshedSurface to shift (1, 0, 0)" << endl;
            Info<< "original" << nl;
            surf.writeStats(Info);
            Info<< endl;

            ModifiableMeshedSurface<face> tsurf(std::move(surf));
            // ModifiableMeshedSurface<face> xtsurf;
            // xtsurf.transfer(surf);

            Info<< "in-progress" << nl;
            surf.writeStats(Info);
            Info<< endl;

            tsurf.storedPoints() += vector(1, 0, 0);

            surf.transfer(tsurf);

            Info<< "updated" << nl;
            surf.writeStats(Info);
            Info<< endl;

            Info<< "modifier" << nl;
            tsurf.writeStats(Info);
            Info<< endl;
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
            surf.writeStats(Info);
            Info<< "Area         : " << sum(surf.magSf()) << nl
                << endl;
        }

        if (optStdout)
        {
            Info<< surf;
        }
        else
        {
            surf.write(exportName);
        }

        if (args.found("surfMesh"))
        {
            Foam::Time runTime
            (
                args.rootPath(),
                args.caseName()
            );

            // start with "constant"
            runTime.setTime(instant(0, runTime.constant()), 0);

            Info<< "runTime.instance() = " << runTime.instance() << endl;
            Info<< "runTime.timeName() = " << runTime.timeName() << endl;

            Info<< "write MeshedSurface 'yetAnother' via proxy as surfMesh"
                << endl;
            surf.write
            (
                runTime,
                "yetAnother"
            );

            surfMesh surfIn
            (
                IOobject
                (
                    "default",
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );


            MeshedSurface<face> surfIn2(runTime, "foobar");

            Info<<"surfIn2 = " << surfIn2.size() << endl;
            Info<< "surfIn = " << surfIn.size() << endl;

            Info<< "writing surfMesh as obj = oldSurfIn.obj" << endl;

            using Foam::surfMesh;
            surfIn.write(fileName("oldSurfIn.obj"));

            Info<< "runTime.instance() = " << runTime.instance() << endl;

            surfMesh surfOut
            (
                IOobject
                (
                    "mySurf",
                    runTime.instance(),
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                std::move(surf)
            );

            Info<< "writing surfMesh as well: " << surfOut.objectPath() << endl;
            surfOut.write();

            surfLabelField zoneIds
            (
                IOobject
                (
                    "zoneIds",
                    surfOut.instance(),
                    surfOut,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surfOut,
                dimless
            );

            Info<<" surf name= " << surfOut.name() <<nl;
            Info<< "rename to anotherSurf" << endl;
            surfOut.rename("anotherSurf");

            Info<<" surf name= " << surfOut.name() <<nl;

            // advance time to 1
            runTime.setTime(instant(1), 1);
            surfOut.setInstance(runTime.timeName());


            Info<< "writing surfMesh again well: " << surfOut.objectPath()
                << endl;
            surfOut.write();

            // write directly
            surfOut.surfMesh::write(fileName("someName.ofs"));

#if 1
            const surfZoneList& zones = surfOut.surfZones();
            forAll(zones, zoneI)
            {
                SubList<label>
                (
                    zoneIds,
                    zones[zoneI].size(),
                    zones[zoneI].start()
                ) = zoneI;
            }

            Info<< "write zoneIds (for testing only): "
                << zoneIds.objectPath() << endl;
            zoneIds.write();

            surfPointLabelField pointIds
            (
                IOobject
                (
                    "zoneIds.",
//                    "pointIds",
                    surfOut.instance(),
//                    "pointFields",
                    surfOut,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surfOut,
                dimless
            );

            forAll(pointIds, i)
            {
                pointIds[i] = i;
            }

            Info<< "write pointIds (for testing only): "
                << pointIds.objectPath() << endl;
            pointIds.write();

            Info<<"surfMesh with these names: " << surfOut.names() << endl;

#endif
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
