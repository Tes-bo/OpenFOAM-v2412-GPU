/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    setsToZones

Group
    grpMeshManipulationUtilities

Description
    Add pointZones/faceZones/cellZones to the mesh from similar named
    pointSets/faceSets/cellSets.

    There is one catch: for faceZones you also need to specify a flip
    condition which basically denotes the side of the face. In this app
    it reads a cellSet (xxxCells if 'xxx' is the name of the faceSet) which
    is the masterCells of the zone.
    There are lots of situations in which this will go wrong but it is the
    best I can think of for now.

    If one is not interested in sideNess specify the -noFlipMap
    command line option.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addNote
    (
        "Add point/face/cell Zones from similarly named point/face/cell Sets"
    );

    argList::addBoolOption
    (
        "noFlipMap",
        "Ignore orientation of faceSet"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const bool noFlipMap = args.found("noFlipMap");

    // Allow override of time from specified time options, or no-op
    timeSelector::setTimeIfPresent(runTime, args);

    #include "createNamedPolyMesh.H"

    const fileName setsSubPath(mesh.meshDir()/"sets");

    // Search for list of objects for the time of the mesh
    word setsInstance = runTime.findInstance
    (
        setsSubPath,
        word::null,
        IOobject::MUST_READ,
        mesh.facesInstance()
    );

    IOobjectList objects(mesh, setsInstance, polyMesh::meshSubDir/"sets");

    Info<< "Searched : " << setsInstance/setsSubPath
        << nl
        << "Found    : " << objects.names() << nl
        << endl;


    for (const IOobject& io : objects.csorted<pointSet>())
    {
        // Not in memory. Load it.
        pointSet set(io);
        labelList pointLabels(set.sortedToc());

        // The original number of zones
        const label nOrigZones = mesh.pointZones().size();

        // Get existing or create new empty zone
        pointZone& zn = mesh.pointZones()(set.name());

        if (nOrigZones == mesh.pointZones().size())
        {
            Info<< "Overwriting contents of existing pointZone "
                << zn.index()
                << " with that of set " << set.name() << "." << endl;
        }
        else
        {
            Info<< "Adding set " << set.name() << " as a pointZone." << endl;
        }

        zn = pointLabels;

        mesh.pointZones().writeOpt(IOobject::AUTO_WRITE);
        mesh.pointZones().instance() = mesh.facesInstance();
    }


    wordHashSet slaveCellSets;

    for (const IOobject& io : objects.csorted<faceSet>())
    {
        // Not in memory. Load it.
        faceSet set(io);
        labelList faceLabels(set.sortedToc());

        DynamicList<label> addressing(set.size());
        DynamicList<bool> flipMap(set.size());

        if (noFlipMap)
        {
            // No flip map.
            forAll(faceLabels, i)
            {
                label facei = faceLabels[i];
                addressing.append(facei);
                flipMap.append(false);
            }
        }
        else
        {
            const word setName(set.name() + "SlaveCells");

            Info<< "Trying to load cellSet " << setName
                << " to find out the slave side of the zone." << nl
                << "If you do not care about the flipMap"
                << " (i.e. do not use the sideness)" << nl
                << "use the -noFlipMap command line option."
                << endl;

            // Load corresponding cells
            cellSet cells(mesh, setName);

            // Store setName to exclude from cellZones further on
            slaveCellSets.insert(setName);

            forAll(faceLabels, i)
            {
                label facei = faceLabels[i];

                bool flip = false;

                if (mesh.isInternalFace(facei))
                {
                    if
                    (
                        cells.found(mesh.faceOwner()[facei])
                    && !cells.found(mesh.faceNeighbour()[facei])
                    )
                    {
                        flip = false;
                    }
                    else if
                    (
                       !cells.found(mesh.faceOwner()[facei])
                     && cells.found(mesh.faceNeighbour()[facei])
                    )
                    {
                        flip = true;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "One of owner or neighbour of internal face "
                            << facei << " should be in cellSet " << cells.name()
                            << " to be able to determine orientation." << endl
                            << "Face:" << facei
                            << " own:" << mesh.faceOwner()[facei]
                            << " OwnInCellSet:"
                            << cells.found(mesh.faceOwner()[facei])
                            << " nei:" << mesh.faceNeighbour()[facei]
                            << " NeiInCellSet:"
                            << cells.found(mesh.faceNeighbour()[facei])
                            << abort(FatalError);
                    }
                }
                else
                {
                    if (cells.found(mesh.faceOwner()[facei]))
                    {
                        flip = false;
                    }
                    else
                    {
                        flip = true;
                    }
                }

                addressing.append(facei);
                flipMap.append(flip);
            }
        }

        // The original number of zones
        const label nOrigZones = mesh.faceZones().size();

        // Get existing or create new empty zone
        faceZone& zn = mesh.faceZones()(set.name());

        if (nOrigZones == mesh.faceZones().size())
        {
            Info<< "Overwriting contents of existing faceZone "
                << zn.index()
                << " with that of set " << set.name() << "." << endl;
        }
        else
        {
            Info<< "Adding set " << set.name() << " as a faceZone." << endl;
        }

        zn.resetAddressing
        (
            addressing.shrink(),
            flipMap.shrink()
        );

        mesh.faceZones().writeOpt(IOobject::AUTO_WRITE);
        mesh.faceZones().instance() = mesh.facesInstance();
    }



    for (const IOobject& io : objects.csorted<cellSet>())
    {
        if (!slaveCellSets.found(io.name()))
        {
            // Not in memory. Load it.
            cellSet set(io);
            labelList cellLabels(set.sortedToc());

            // The original number of zones
            const label nOrigZones = mesh.cellZones().size();

            cellZone& zn = mesh.cellZones()(set.name());

            if (nOrigZones == mesh.cellZones().size())
            {
                Info<< "Overwriting contents of existing cellZone "
                    << zn.index()
                    << " with that of set " << set.name() << "." << endl;
            }
            else
            {
                Info<< "Adding set " << set.name() << " as a cellZone." << endl;
            }

            zn = cellLabels;

            mesh.cellZones().writeOpt(IOobject::AUTO_WRITE);
            mesh.cellZones().instance() = mesh.facesInstance();
        }
    }


    Info<< "Writing mesh." << endl;

    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
