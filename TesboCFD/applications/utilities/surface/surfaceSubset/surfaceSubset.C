/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    surfaceSubset

Group
    grpSurfaceUtilities

Description
    A surface analysis tool that subsets the triSurface to choose a
    region of interest. Based on subsetMesh.

\*---------------------------------------------------------------------------*/

#include "triSurfaceSearch.H"
#include "MeshedSurfaces.H"
#include "argList.H"
#include "Fstream.H"
#include "IOdictionary.H"
#include "boundBox.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "Random.H"
#include "volumeType.H"
#include "plane.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A surface analysis tool that subsets the surface to choose a"
        " region of interest."
    );

    argList::noParallel();
    argList::addArgument("dict", "The surfaceSubsetDict");
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");
    argList args(argc, argv);

    Info<< "Reading dictionary " << args[1] << " ..." << endl;
    IFstream dictFile(args.get<fileName>(1));
    dictionary meshSubsetDict(dictFile);

    Info<< "Reading surface " << args[2] << " ..." << endl;
    meshedSurface surf1(args.get<fileName>(2));

    const auto outFileName(args.get<fileName>(3));

    Info<< "Original:" << endl;
    surf1.writeStats(Info);
    Info<< endl;


    labelList markedPoints
    (
        meshSubsetDict.lookup("localPoints")
    );

    labelList markedEdges
    (
        meshSubsetDict.lookup("edges")
    );

    labelList markedFaces
    (
        meshSubsetDict.lookup("faces")
    );

    pointField markedZone
    (
        meshSubsetDict.lookup("zone")
    );


    boundBox zoneBb;

    if (markedZone.size())
    {
        if (markedZone.size() != 2)
        {
            FatalErrorInFunction
                << "zone specification should be two points, min and max of "
                << "the boundingbox" << endl
                << "zone:" << markedZone
                << exit(FatalError);
        }

        zoneBb.min() = markedZone[0];
        zoneBb.max() = markedZone[1];

        if (!zoneBb.good())
        {
            WarningInFunction
                << "Defined zone is invalid: " << zoneBb << nl;
        }
    }


    const bool addFaceNeighbours =
        meshSubsetDict.get<bool>("addFaceNeighbours");

    const bool invertSelection =
        meshSubsetDict.getOrDefault("invertSelection", false);

    // Mark the cells for the subset

    // Faces to subset
    bitSet facesToSubset(surf1.size(), false);


    //
    // Faces connected to "localPoints"
    //

    if (markedPoints.size())
    {
        Info<< "Found " << markedPoints.size() << " marked point(s)." << endl;

        for (const label pointi : markedPoints)
        {
            if (pointi < 0 || pointi >= surf1.nPoints())
            {
                FatalErrorInFunction
                    << "localPoint label " << pointi << "out of range."
                    << " Surface has " << surf1.nPoints() << " localPoints."
                    << exit(FatalError);
            }

            const labelList& curFaces = surf1.pointFaces()[pointi];

            facesToSubset.set(curFaces);
        }
    }


    //
    // Faces connected to "edges"
    //

    if (markedEdges.size())
    {
        Info<< "Found " << markedEdges.size() << " marked edge(s)." << endl;

        for (const label edgei : markedEdges)
        {
            if (edgei < 0 || edgei >= surf1.nEdges())
            {
                FatalErrorInFunction
                    << "edge label " << edgei << "out of range."
                    << " Surface has " << surf1.nEdges() << " edges."
                    << exit(FatalError);
            }

            const labelList& curFaces = surf1.edgeFaces()[edgei];

            facesToSubset.set(curFaces);
        }
    }


    //
    // Faces with centre inside "zone"
    //

    if (zoneBb.good())
    {
        Info<< "Using zone " << zoneBb << endl;

        forAll(surf1, facei)
        {
            const point centre = surf1[facei].centre(surf1.points());

            if (zoneBb.contains(centre))
            {
                facesToSubset.set(facei);
            }
        }
    }


    //
    // Faces on certain side of surface
    //

    if (meshSubsetDict.found("surface"))
    {
        const dictionary& surfDict = meshSubsetDict.subDict("surface");

        const auto surfName(surfDict.get<fileName>("name"));

        const volumeType::type volType =
        (
            surfDict.getOrDefault("outside", false)
          ? volumeType::OUTSIDE
          : volumeType::INSIDE
        );

        Info<< "Selecting faces with centre located "
            << volumeType::names[volType] << " of surface "
            << surfName << endl;

        // Read surface to select on
        triSurface selectSurf(surfName);

        triSurfaceSearch searchSelectSurf
        (
            selectSurf,
            indexedOctree<treeDataTriSurface>::perturbTol(),
            8
        );

        const indexedOctree<treeDataTriSurface>& selectTree =
            searchSelectSurf.tree();

        // Check if face (centre) is in outside or inside.
        forAll(surf1, facei)
        {
            if (!facesToSubset[facei])
            {
                const point fc(surf1[facei].centre(surf1.points()));

                if (volType == selectTree.getVolumeType(fc))
                {
                    facesToSubset.set(facei);
                }
            }
        }
    }


    if (meshSubsetDict.found("plane"))
    {
        const dictionary& planeDict = meshSubsetDict.subDict("plane");

        const plane pl(planeDict);
        const scalar distance(planeDict.get<scalar>("distance"));
        const scalar cosAngle(planeDict.get<scalar>("cosAngle"));

        // Select all triangles that are close to the plane and
        // whose normal aligns with the plane as well.

        forAll(surf1.faceCentres(), facei)
        {
            const point& fc = surf1.faceCentres()[facei];
            const point& nf = surf1.faceNormals()[facei];

            if (pl.distance(fc) < distance && mag(pl.normal() & nf) > cosAngle)
            {
                facesToSubset.set(facei);
            }
        }
    }



    //
    // pick up specified "faces"
    //

    // Number of additional faces picked up because of addFaceNeighbours
    label nFaceNeighbours = 0;

    if (markedFaces.size())
    {
        Info<< "Found " << markedFaces.size() << " marked face(s)." << endl;

        // Check and mark faces to pick up
        for (const label facei : markedFaces)
        {
            if (facei < 0 || facei >= surf1.size())
            {
                FatalErrorInFunction
                    << "Face label " << facei << "out of range."
                    << " Surface has " << surf1.size() << " faces."
                    << exit(FatalError);
            }

            // Mark the face
            facesToSubset.set(facei);

            // Mark its neighbours if requested
            if (addFaceNeighbours)
            {
                const labelList& curFaces = surf1.faceFaces()[facei];

                for (const label neiFacei : curFaces)
                {
                    if (facesToSubset.set(neiFacei))
                    {
                        ++nFaceNeighbours;
                    }
                }
            }
        }
    }

    if (addFaceNeighbours)
    {
        Info<< "Added " << nFaceNeighbours
            << " faces because of addFaceNeighbours" << endl;
    }


    if (invertSelection)
    {
        Info<< "Inverting selection." << endl;

        facesToSubset.flip();
    }


    // Create subsetted surface
    meshedSurface surf2(surf1.subsetMesh(facesToSubset));

    Info<< "Subset:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    Info<< "Writing surface to " << outFileName << endl;

    surf2.write(outFileName);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
