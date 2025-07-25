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
    surfaceAdd

Group
    grpSurfaceUtilities

Description
    Add two surfaces. Does geometric merge on points.
    Does not check for overlapping/intersecting triangles.

    Keeps patches separate by renumbering.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "Fstream.H"
#include "triFace.H"
#include "triFaceList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Add two surfaces via a geometric merge on points."
        " Does not check for overlapping/intersecting triangles."
    );

    argList::noParallel();
    argList::addArgument("surface1", "The input surface file 1");
    argList::addArgument("surface2", "The input surface file 2");
    argList::addArgument("output", "The output surface file");

    argList::addOption
    (
        "points",
        "file",
        "Provide additional points"
    );
    argList::addBoolOption
    (
        "mergeRegions",
        "Combine regions from both surfaces"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor on input surfaces"
    );
    argList::addVerboseOption();

    argList args(argc, argv);

    const auto inFileName1 = args.get<fileName>(1);
    const auto inFileName2 = args.get<fileName>(2);
    const auto outFileName = args.get<fileName>(3);

    const int  optVerbose   = args.verbose();
    const bool addPoint     = args.found("points");
    const bool mergeRegions = args.found("mergeRegions");

    const scalar scaleFactor = args.getOrDefault<scalar>("scale", -1);

    if (addPoint)
    {
        Info<< "Reading a surface and adding points from a file"
            << "; merging the points and writing the surface to another file"
            << nl << endl;

        Info<< "Surface  : " << inFileName1<< nl
            << "Points   : " << args.get<fileName>("points") << nl
            << "Writing  : " << outFileName << nl << endl;
    }
    else
    {
        Info<< "Reading two surfaces"
            << "; merging points and writing the surface to another file"
            << nl << endl;

        if (mergeRegions)
        {
            Info<< "Regions from the two files will get merged" << nl
                << "Do not use this option if you want to keep the regions"
                << " separate" << nl << endl;
        }
        else
        {
            Info<< "Regions from the two files will not get merged" << nl
                << "Regions from " << inFileName2 << " will get offset so"
                << " as not to overlap with the regions in " << inFileName1
                << nl << endl;
        }


        Info<< "Surface1 : " << inFileName1<< nl
            << "Surface2 : " << inFileName2<< nl
            << "Writing  : " << outFileName << nl << endl;
    }

    if (scaleFactor > 0)
    {
        Info<< "Scaling  : " << scaleFactor << nl;
    }

    const triSurface surface1(inFileName1, scaleFactor);
    Info<< "Surface1:" << endl;
    surface1.writeStats(Info);
    Info<< endl;

    const pointField& points1 = surface1.points();

    // Final surface
    triSurface combinedSurf;

    if (addPoint)
    {
        IFstream pointsFile(args.get<fileName>("points"));
        const pointField extraPoints(pointsFile);

        Info<< "Additional Points:" << extraPoints.size() << endl;

        vectorField pointsAll(points1);
        label pointi = pointsAll.size();
        pointsAll.setSize(pointsAll.size() + extraPoints.size());

        for (const auto& pt : extraPoints)
        {
            pointsAll[pointi++] = pt;
        }

        combinedSurf = triSurface(surface1, surface1.patches(), pointsAll);
    }
    else
    {
        const triSurface surface2(inFileName2, scaleFactor);
        Info<< "Surface2:" << endl;
        surface2.writeStats(Info);
        Info<< endl;


        // Make new storage
        List<labelledTri> facesAll(surface1.size() + surface2.size());

        const pointField& points2 = surface2.points();

        vectorField pointsAll(points1.size() + points2.size());


        label pointi = 0;
        // Copy points1 into pointsAll
        for (const auto& pt : points1)
        {
            pointsAll[pointi++] = pt;
        }
        // Add surface2 points
        for (const auto& pt : points2)
        {
            pointsAll[pointi++] = pt;
        }


        label trianglei = 0;

        // Determine map for both regions
        label nNewPatches = 0;
        labelList patch1Map(surface1.patches().size());
        labelList patch2Map(surface2.patches().size());

        if (mergeRegions)
        {
            HashTable<label> nameToPatch;

            forAll(surface1.patches(), i)
            {
                const word& name = surface1.patches()[i].name();

                // Lookup or insert
                const label combinedi = nameToPatch(name, nameToPatch.size());

                patch1Map[i] = combinedi;
            }

            // Determine map for surface2 regions

            forAll(surface2.patches(), i)
            {
                const word& name = surface2.patches()[i].name();

                // Lookup or insert
                const label combinedi = nameToPatch(name, nameToPatch.size());

                patch2Map[i] = combinedi;
            }

            nNewPatches = nameToPatch.size();
        }
        else
        {
            Info<< "Surface " << inFileName1
                << " has " << surface1.patches().size()
                << " regions"
                << nl
                << "All region numbers in " << inFileName2 << " will be offset"
                << " by this amount" << nl << endl;

            patch1Map = identity(surface1.patches().size());
            patch2Map = identity(surface2.patches().size(), patch1Map.size());

            nNewPatches = surface1.patches().size()+surface2.patches().size();
        }


        // Copy triangles1 into trianglesAll
        for (const labelledTri& tri : surface1)
        {
            labelledTri& destTri = facesAll[trianglei++];

            destTri.triFace::operator=(tri);
            destTri.region() = patch1Map[tri.region()];
        }

        // Add (renumbered) surface2 triangles
        for (const labelledTri& tri : surface2)
        {
            labelledTri& destTri = facesAll[trianglei++];
            destTri[0] = tri[0] + points1.size();
            destTri[1] = tri[1] + points1.size();
            destTri[2] = tri[2] + points1.size();
            destTri.region() = patch2Map[tri.region()];
        }


        geometricSurfacePatchList newPatches(nNewPatches);
        forAll(surface1.patches(), patchi)
        {
            newPatches[patch1Map[patchi]] = surface1.patches()[patchi];
        }
        forAll(surface2.patches(), patchi)
        {
            newPatches[patch2Map[patchi]] = surface2.patches()[patchi];
        }

        Info<< "New patches:" << nl;
        forAll(newPatches, patchi)
        {
            Info<< "    " << patchi << '\t' << newPatches[patchi].name() << nl;
        }
        Info<< endl;


        // Construct new surface mesh
        combinedSurf = triSurface(facesAll, newPatches, pointsAll);
    }

    // Merge all common points and do some checks
    combinedSurf.cleanup(optVerbose);

    Info<< "Merged surface:" << endl;

    combinedSurf.writeStats(Info);

    Info<< endl;

    Info<< "Writing : " << outFileName << endl;

    // If merging regions also sort
    combinedSurf.write(outFileName, mergeRegions);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
