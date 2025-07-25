/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    topoSet

Group
    grpMeshManipulationUtilities

Description
    Operates on cellSets/faceSets/pointSets through a dictionary,
    normally system/topoSetDict

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "topoSetSource.H"
#include "globalMeshData.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "pointZoneSet.H"
#include "IOdictionary.H"
#include "namedDictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void printMesh(const Time& runTime, const polyMesh& mesh)
{
    Info<< "Time:" << runTime.timeName()
        << "  cells:" << mesh.globalData().nTotalCells()
        << "  faces:" << mesh.globalData().nTotalFaces()
        << "  points:" << mesh.globalData().nTotalPoints()
        << "  patches:" << mesh.boundaryMesh().size()
        << "  bb:" << mesh.bounds() << nl;
}


template<class ZoneType>
void removeZone
(
    ZoneMesh<ZoneType, polyMesh>& zones,
    const word& setName
)
{
    label zoneID = zones.findZoneID(setName);

    if (zoneID != -1)
    {
        Info<< "Removing zone " << setName << " at index " << zoneID << endl;
        // Shuffle to last position
        labelList oldToNew(zones.size());
        label newI = 0;
        forAll(oldToNew, i)
        {
            if (i != zoneID)
            {
                oldToNew[i] = newI++;
            }
        }
        oldToNew[zoneID] = newI;
        zones.reorder(oldToNew);
        // Remove last element
        zones.setSize(zones.size()-1);
        zones.clearAddressing();
        if (!zones.write())
        {
            WarningInFunction << "Failed writing zone " << setName << endl;
        }
        fileHandler().flush();
    }
}


// Physically remove a set
void removeSet
(
    const polyMesh& mesh,
    const word& setType,
    const word& setName
)
{
    // Remove the file
    IOobjectList objects
    (
        mesh,
        mesh.time().findInstance
        (
            polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::READ_IF_PRESENT,
            mesh.facesInstance()
        ),
        polyMesh::meshSubDir/"sets"
    );

    if (objects.found(setName))
    {
        // Remove file
        fileName object = objects[setName]->objectPath();
        Info<< "Removing file " << object << endl;
        rm(object);
    }

    // See if zone
    if (setType == cellZoneSet::typeName)
    {
        removeZone
        (
            const_cast<cellZoneMesh&>(mesh.cellZones()),
            setName
        );
    }
    else if (setType == faceZoneSet::typeName)
    {
        removeZone
        (
            const_cast<faceZoneMesh&>(mesh.faceZones()),
            setName
        );
    }
    else if (setType == pointZoneSet::typeName)
    {
        removeZone
        (
            const_cast<pointZoneMesh&>(mesh.pointZones()),
            setName
        );
    }
}


polyMesh::readUpdateState meshReadUpdate(polyMesh& mesh)
{
    polyMesh::readUpdateState stat = mesh.readUpdate();

    switch(stat)
    {
        case polyMesh::UNCHANGED:
        {
            Info<< "    mesh not changed." << endl;
            break;
        }
        case polyMesh::POINTS_MOVED:
        {
            Info<< "    points moved; topology unchanged." << endl;
            break;
        }
        case polyMesh::TOPO_CHANGE:
        {
            Info<< "    topology changed; patches unchanged." << nl
                << "    ";
            printMesh(mesh.time(), mesh);
            break;
        }
        case polyMesh::TOPO_PATCH_CHANGE:
        {
            Info<< "    topology changed and patches changed." << nl
                << "    ";
            printMesh(mesh.time(), mesh);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Illegal mesh update state "
                << stat  << abort(FatalError);
            break;
        }
    }
    return stat;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Operates on cellSets/faceSets/pointSets through a dictionary,"
        " normally system/topoSetDict"
    );

    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    argList::addOption("dict", "file", "Alternative topoSetDict");

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "noSync",
        "Do not synchronise selection across coupled patches"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedPolyMesh.H"

    const bool noSync = args.found("noSync");

    const word dictName("topoSetDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary topoSetDict(dictIO);

    // Read set construct info from dictionary
    List<namedDictionary> actionEntries(topoSetDict.lookup("actions"));

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        // Optionally re-read mesh
        meshReadUpdate(mesh);

        // Execute all actions
        for (const namedDictionary& actionEntry : actionEntries)
        {
            const dictionary& dict = actionEntry.dict();
            if (dict.empty())
            {
                continue;
            }
            const word setName(dict.get<word>("name"));
            const word setType(dict.get<word>("type"));

            const topoSetSource::setAction action =
                topoSetSource::actionNames.get("action", dict);

            autoPtr<topoSet> currentSet;

            switch (action)
            {
                case topoSetSource::NEW :
                case topoSetSource::CLEAR :
                {
                    currentSet = topoSet::New(setType, mesh, setName, 16384);
                    Info<< "Created "
                        << currentSet().type() << ' ' << setName << endl;
                    break;
                }

                case topoSetSource::IGNORE :
                    continue;  // Nothing to do
                    break;

                case topoSetSource::REMOVE :
                    // Nothing to load
                    break;

                default:
                {
                    // Load set
                    currentSet = topoSet::New
                    (
                        setType,
                        mesh,
                        setName,
                        IOobject::MUST_READ
                    );

                    Info<< "Read set "
                        << currentSet().type() << ' ' << setName
                        << " size:"
                        << returnReduce(currentSet().size(), sumOp<label>())
                        << endl;
                }
            }

            // Handle special actions (clear, invert) locally,
            // the other actions through sources.
            switch (action)
            {
                case topoSetSource::NEW :
                case topoSetSource::ADD :
                case topoSetSource::SUBTRACT :
                {
                    const word sourceType(dict.get<word>("source"));

                    Info<< "    Applying source " << sourceType << endl;
                    autoPtr<topoSetSource> source = topoSetSource::New
                    (
                        sourceType,
                        mesh,
                        dict.optionalSubDict("sourceInfo")
                    );

                    source().applyToSet(action, currentSet());
                    // Synchronize for coupled patches.
                    if (!noSync) currentSet().sync(mesh);
                    if (!currentSet().write())
                    {
                        WarningInFunction
                            << "Failed writing set "
                            << currentSet().objectPath() << endl;
                    }
                    fileHandler().flush();
                    break;
                }

                case topoSetSource::SUBSET :
                {
                    const word sourceType(dict.get<word>("source"));

                    Info<< "    Applying source " << sourceType << endl;
                    autoPtr<topoSetSource> source = topoSetSource::New
                    (
                        sourceType,
                        mesh,
                        dict.optionalSubDict("sourceInfo")
                    );

                    // Backup current set.
                    autoPtr<topoSet> oldSet
                    (
                        topoSet::New
                        (
                            setType,
                            mesh,
                            currentSet().name() + "_old2",
                            currentSet()
                        )
                    );

                    currentSet().clear();
                    source().applyToSet(topoSetSource::NEW, currentSet());

                    // Combine new value of currentSet with old one.
                    currentSet().subset(oldSet());
                    // Synchronize for coupled patches.
                    if (!noSync) currentSet().sync(mesh);
                    if (!currentSet().write())
                    {
                        WarningInFunction
                            << "Failed writing set "
                            << currentSet().objectPath() << endl;
                    }
                    fileHandler().flush();

                    break;
                }

                case topoSetSource::CLEAR :
                {
                    Info<< "    Clearing " << currentSet().type() << endl;
                    currentSet().clear();
                    if (!currentSet().write())
                    {
                        WarningInFunction
                            << "Failed writing set "
                            << currentSet().objectPath() << endl;
                    }
                    fileHandler().flush();

                    break;
                }

                case topoSetSource::INVERT :
                {
                    Info<< "    Inverting " << currentSet().type() << endl;
                    currentSet().invert(currentSet().maxSize(mesh));
                    if (!currentSet().write())
                    {
                        WarningInFunction
                            << "Failed writing set "
                            << currentSet().objectPath() << endl;
                    }
                    fileHandler().flush();

                    break;
                }

                case topoSetSource::REMOVE :
                {
                    Info<< "    Removing set" << endl;
                    removeSet(mesh, setType, setName);

                    break;
                }

                default:
                    WarningInFunction
                        << "Unhandled action: "
                        << topoSetSource::actionNames[action] << endl;
            }

            if (currentSet)
            {
                Info<< "    "
                    << currentSet().type() << ' '
                    << currentSet().name() << " now size "
                    << returnReduce(currentSet().size(), sumOp<label>())
                    << endl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
