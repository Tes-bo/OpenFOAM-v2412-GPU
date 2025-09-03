/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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
    setSet

Group
    grpMeshManipulationUtilities

Description
    Manipulate a cell/face/point Set or Zone interactively.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "StringStream.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "topoSetSource.H"
#include "Fstream.H"
#include "foamVtkWriteTopoSet.H"
#include "IOobjectList.H"
#include "cellZoneSet.H"
#include "faceZoneSet.H"
#include "pointZoneSet.H"
#include "timeSelector.H"

#include <stdio.h>

#ifdef HAVE_LIBREADLINE
    #include <readline/readline.h>
    #include <readline/history.h>

    static const char* historyFile = ".setSet";
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Write set to VTK readable files
void writeVTK
(
    const polyMesh& mesh,
    const topoSet& currSet,
    const fileName& outputName
)
{
    if
    (
        !vtk::writeTopoSet
        (
            mesh,
            currSet,
            vtk::formatType::INLINE_BASE64,     // XML-binary
            // vtk::formatType::LEGACY_BINARY,
            outputName,
            false // Not parallel
        )
    )
    {
        WarningInFunction
            << "Don't know how to handle set of type "
            << currSet.type() << nl;
    }
}


void printHelp(Ostream& os)
{
    os  << "Please type 'help', 'list', 'quit', 'time ddd'"
        << " or a set command after prompt." << nl
        << "'list' will show all current cell/face/point sets." << nl
        << "'time ddd' will change the current time." << nl
        << nl
        << "A set command should be of the following form" << nl
        << nl
        << "    cellSet|faceSet|pointSet <setName> <action> <source>"
        << nl
        << nl
        << "The <action> is one of" << nl
        << "    list            - prints the contents of the set" << nl
        << "    clear           - clears the set" << nl
        << "    invert          - inverts the set" << nl
        << "    remove          - remove the set" << nl
        << "    new <source>    - use all elements from the source set" << nl
        << "    add <source>    - adds all elements from the source set" << nl
        << "    subtract <source> - subtract the source set elements" << nl
        << "    subset <source> - combines current set with the source set"
        << nl
        << nl
        << "The sources come in various forms. Type a wrong source"
        << " to see all the types available." << nl
        << nl
        << "Example: pick up all cells connected by point or face to patch"
        << " movingWall" << nl
        << nl
        << "Pick up all faces of patch:" << nl
        << "    faceSet f0 new patchToFace movingWall" << nl
        << "Add faces 0,1,2:" << nl
        << "    faceSet f0 add labelToFace (0 1 2)" << nl
        << "Pick up all points used by faces in faceSet f0:" << nl
        << "    pointSet p0 new faceToPoint f0 all" << nl
        << "Pick up cell which has any face in f0:" << nl
        << "    cellSet c0 new faceToCell f0 any" << nl
        << "Add cells which have any point in p0:" << nl
        << "    cellSet c0 add pointToCell p0 any" << nl
        << "List set:" << nl
        << "    cellSet c0 list" << nl
        << nl
        << "Zones can be set using zoneSets from corresponding sets:" << nl
        << "    cellZoneSet c0Zone new setToCellZone c0" << nl
        << "    faceZoneSet f0Zone new setToFaceZone f0" << nl
        << nl
        << "or if orientation is important:" << nl
        << "    faceZoneSet f0Zone new setsToFaceZone f0 c0" << nl
        << nl
        << "ZoneSets can be manipulated using the general actions:" << nl
        << "    list            - prints the contents of the set" << nl
        << "    clear           - clears the set" << nl
        << "    invert          - inverts the set (undefined orientation)"
        << nl
        << "    remove          - remove the set" << nl
        << endl;
}


template<class SetType>
void printSets(Ostream& os, const IOobjectList& objects)
{
    label n = 0;

    for (const IOobject& io : objects.csorted<SetType>())
    {
        SetType set(io);
        if (!n++) os << SetType::typeName << "s:" << nl;
        os  << '\t' << set.name() << "\tsize:" << set.size() << endl;
    }
}


template<class ZoneType>
void printZones(Ostream& os, const ZoneMesh<ZoneType, polyMesh>& zones)
{
    label n = 0;

    for (const ZoneType& zn : zones)
    {
        if (!n++) os << ZoneType::typeName << "s:" << nl;
        os  << '\t' << zn.name() << "\tsize:" << zn.size() << endl;
    }
}


void printAllSets(const polyMesh& mesh, Ostream& os)
{
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

    printSets<cellSet>(os, objects);
    printSets<faceSet>(os, objects);
    printSets<pointSet>(os, objects);

    printZones(os, mesh.cellZones());
    printZones(os, mesh.faceZones());
    printZones(os, mesh.pointZones());

    os  << endl;
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
        zones.write();
        // Force flushing so we know it has finished writing
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


// Read command and execute. Return true if ok, false otherwise.
bool doCommand
(
    const polyMesh& mesh,
    const word& setType,
    const word& setName,
    const word& actionName,
    const bool writeVTKFile,
    const bool writeCurrentTime,
    const bool noSync,
    Istream& is
)
{
    // Get some size estimate for set.
    const globalMeshData& parData = mesh.globalData();

    label typSize =
        max
        (
            parData.nTotalCells(),
            max
            (
                parData.nTotalFaces(),
                parData.nTotalPoints()
            )
        )
      / (10*Pstream::nProcs());


    bool ok = true;

    // Set to work on
    autoPtr<topoSet> currentSetPtr;

    word sourceType;

    try
    {
        topoSetSource::setAction action =
            topoSetSource::actionNames[actionName];

        switch (action)
        {
            case topoSetSource::REMOVE :
            {
                removeSet(mesh, setType, setName);
                break;
            }

            case topoSetSource::NEW :
            case topoSetSource::CLEAR :
            {
                currentSetPtr = topoSet::New(setType, mesh, setName, typSize);
                break;
            }

            case topoSetSource::IGNORE :
                // Nothing to do
                break;

            default:
            {
                currentSetPtr = topoSet::New
                (
                    setType,
                    mesh,
                    setName,
                    IOobject::MUST_READ
                );

                topoSet& currentSet = currentSetPtr();
                // Presize it according to current mesh data.
                currentSet.reserve(max(currentSet.size(), typSize));
            }
        }

        if (currentSetPtr)
        {
            topoSet& currentSet = *currentSetPtr;

            Info<< "    Set:" << currentSet.name()
                << "  Size:" << returnReduce(currentSet.size(), sumOp<label>())
                << "  Action:" << actionName
                << endl;

            switch (action)
            {
                case topoSetSource::CLEAR :
                {
                    // Already handled above by not reading
                    break;
                }

                case topoSetSource::INVERT :
                {
                    currentSet.invert(currentSet.maxSize(mesh));
                    break;
                }

                case topoSetSource::LIST :
                {
                    currentSet.writeDebug(Pout, mesh, 100);
                    Pout<< endl;
                    break;
                }

                case topoSetSource::SUBSET :
                {
                    if (is >> sourceType)
                    {
                        autoPtr<topoSetSource> setSource
                        (
                            topoSetSource::New
                            (
                                sourceType,
                                mesh,
                                is
                            )
                        );

                        // Backup current set.
                        autoPtr<topoSet> oldSet
                        (
                            topoSet::New
                            (
                                setType,
                                mesh,
                                currentSet.name() + "_old2",
                                currentSet
                            )
                        );

                        currentSet.clear();
                        setSource().applyToSet(topoSetSource::NEW, currentSet);

                        // Combine new value of currentSet with old one.
                        currentSet.subset(oldSet());
                    }
                    break;
                }

                default:
                {
                    if (is >> sourceType)
                    {
                        autoPtr<topoSetSource> setSource
                        (
                            topoSetSource::New
                            (
                                sourceType,
                                mesh,
                                is
                            )
                        );

                        setSource().applyToSet(action, currentSet);
                    }
                }
            }

            if (action != topoSetSource::LIST)
            {
                // Set will have been modified.

                // Synchronize for coupled patches.
                if (!noSync) currentSet.sync(mesh);

                // Write
                Info<< "    Writing " << currentSet.name()
                    << " (size "
                    << returnReduce(currentSet.size(), sumOp<label>())
                    << ") to "
                    << (
                          currentSet.instance()/currentSet.local()
                        / currentSet.name()
                       );


                if (writeVTKFile)
                {
                    fileName outputName
                    (
                        mesh.time().path()/"VTK"/currentSet.name()
                      / currentSet.name() + "_"
                      + Foam::name(mesh.time().timeIndex())
                    );
                    mkDir(outputName.path());

                    Info<< " and to vtk file "
                        << outputName.relative(mesh.time().path())
                        << nl << nl;

                    writeVTK(mesh, currentSet, outputName);
                }
                else
                {
                    Info<< nl << nl;
                }

                if (writeCurrentTime)
                {
                    currentSet.instance() = mesh.time().timeName();
                }
                if (!currentSet.write())
                {
                    WarningInFunction
                        << "Failed writing set "
                        << currentSet.objectPath() << endl;
                }
                // Make sure writing is finished
                fileHandler().flush();
            }
        }
    }
    catch (const Foam::IOerror& fIOErr)
    {
        ok = false;

        Pout<< fIOErr.message().c_str() << endl;

        if (sourceType.size())
        {
            Pout<< topoSetSource::usage(sourceType).c_str();
        }
    }
    catch (const Foam::error& fErr)
    {
        ok = false;

        Pout<< fErr.message().c_str() << endl;

        if (sourceType.size())
        {
            Pout<< topoSetSource::usage(sourceType).c_str();
        }
    }

    return ok;
}


// Status returned from parsing the first token of the line
enum commandStatus
{
    QUIT,           // quit program
    INVALID,        // token is not a valid set manipulation command
    VALIDSETCMD,    // ,,    is a valid     ,,
    VALIDZONECMD    // ,,    is a valid     zone      ,,
};


void printMesh(const Time& runTime, const polyMesh& mesh)
{
    Info<< "Time:" << runTime.timeName()
        << "  cells:" << mesh.globalData().nTotalCells()
        << "  faces:" << mesh.globalData().nTotalFaces()
        << "  points:" << mesh.globalData().nTotalPoints()
        << "  patches:" << mesh.boundaryMesh().size()
        << "  bb:" << mesh.bounds() << nl;
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


commandStatus parseType
(
    Time& runTime,
    polyMesh& mesh,
    const word& setType,
    IStringStream& is
)
{
    if (setType.empty())
    {
        Info<< "Type 'help' for usage information" << endl;

        return INVALID;
    }
    else if (setType == "help")
    {
        printHelp(Info);

        return INVALID;
    }
    else if (setType == "list")
    {
        printAllSets(mesh, Info);

        return INVALID;
    }
    else if (setType == "time")
    {
        scalar requestedTime = readScalar(is);
        instantList Times = runTime.times();

        label nearestIndex = Time::findClosestTimeIndex(Times, requestedTime);

        Info<< "Changing time from " << runTime.timeName()
            << " to " << Times[nearestIndex].name()
            << endl;

        // Set time
        runTime.setTime(Times[nearestIndex], nearestIndex);
        // Optionally re-read mesh
        meshReadUpdate(mesh);

        return INVALID;
    }
    else if (setType == "quit")
    {
        Info<< "Quitting ..." << endl;

        return QUIT;
    }
    else if
    (
        setType == "cellSet"
     || setType == "faceSet"
     || setType == "pointSet"
    )
    {
        return VALIDSETCMD;
    }
    else if
    (
        setType == "cellZoneSet"
     || setType == "faceZoneSet"
     || setType == "pointZoneSet"
    )
    {
        return VALIDZONECMD;
    }
    else
    {
        SeriousErrorInFunction
            << "Illegal command " << setType << endl
            << "Should be one of 'help', 'list', 'time' or a set type :"
            << " 'cellSet', 'faceSet', 'pointSet', 'faceZoneSet'"
            << endl;

        return INVALID;
    }
}


commandStatus parseAction(const word& actionName)
{
    return
    (
        actionName.size() && topoSetSource::actionNames.found(actionName)
      ? VALIDSETCMD : INVALID
    );
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Manipulate a cell/face/point Set or Zone interactively."
    );

    // Specific to topoSet/setSet: quite often we want to block upon writing
    // a set so we can immediately re-read it. So avoid use of threading
    // for set writing.

    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    #include "addRegionOption.H"
    argList::addBoolOption("noVTK", "Do not write VTK files");
    argList::addBoolOption("loop", "Execute batch commands for all timesteps");
    argList::addOption
    (
        "batch",
        "file",
        "Process in batch mode, using input from specified file"
    );
    argList::addBoolOption
    (
        "noSync",
        "Do not synchronise selection across coupled patches"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool writeVTK = !args.found("noVTK");
    const bool loop = args.found("loop");
    const bool batch = args.found("batch");
    const bool noSync = args.found("noSync");

    if (loop && !batch)
    {
        FatalErrorInFunction
            << "Can only loop in batch mode."
            << exit(FatalError);
    }


    #include "createNamedPolyMesh.H"

    // Print some mesh info
    printMesh(runTime, mesh);

    // Print current sets
    printAllSets(mesh, Info);

    // Read history if interactive
    #ifdef HAVE_LIBREADLINE
    if (!batch && !read_history((runTime.path()/historyFile).c_str()))
    {
        Info<< "Successfully read history from " << historyFile << endl;
    }
    #endif


    // Exit status
    int status = 0;


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Handle geometry/topology changes
        meshReadUpdate(mesh);


        // Main command read & execute loop
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<IFstream> fileStreamPtr;

        if (batch)
        {
            const auto batchFile = args.get<fileName>("batch");

            Info<< "Reading commands from file " << batchFile << endl;

            // we cannot handle .gz files
            if (!isFile(batchFile, false))
            {
                FatalErrorInFunction
                    << "Cannot open file " << batchFile << exit(FatalError);
            }

            fileStreamPtr.reset(new IFstream(batchFile));
        }

        Info<< "Please type 'help', 'quit' or a set command after prompt."
            << endl;

        // Whether to quit
        bool quit = false;

        FatalError.throwing(true);
        FatalIOError.throwing(true);

        do
        {
            string rawLine;

            // Type: cellSet, faceSet, pointSet
            word setType;
            // Name of destination set.
            word setName;
            // Action (new, invert etc.)
            word actionName;

            commandStatus stat = INVALID;

            if (fileStreamPtr)
            {
                if (!fileStreamPtr->good())
                {
                    Info<< "End of batch file" << endl;
                    // No error.
                    break;
                }

                fileStreamPtr().getLine(rawLine);

                if (rawLine.size())
                {
                    Info<< "Doing:" << rawLine << endl;
                }
            }
            else
            {
                #ifdef HAVE_LIBREADLINE
                {
                    char* linePtr = readline("readline>");

                    if (linePtr)
                    {
                        rawLine = string(linePtr);

                        if (*linePtr)
                        {
                            add_history(linePtr);
                            write_history(historyFile);
                        }

                        free(linePtr);   // readline uses malloc, not new.
                    }
                    else
                    {
                        break;
                    }
                }
                #else
                {
                    if (!std::cin.good())
                    {
                        Info<< "End of cin" << endl;
                        // No error.
                        break;
                    }
                    Info<< "Command>" << flush;
                    std::getline(std::cin, rawLine);
                }
                #endif
            }

            // Strip off anything after #
            string::size_type i = rawLine.find('#');
            if (i != string::npos)
            {
                rawLine.resize(i);
            }

            if (rawLine.empty())
            {
                continue;
            }

            IStringStream is(rawLine + ' ');

            // Type: cellSet, faceSet, pointSet, faceZoneSet
            is  >> setType;

            stat = parseType(runTime, mesh, setType, is);

            if (stat == VALIDSETCMD || stat == VALIDZONECMD)
            {
                if (is >> setName)
                {
                    if (is >> actionName)
                    {
                        stat = parseAction(actionName);
                    }
                }
            }

            if (stat == QUIT)
            {
                // Make sure to quit
                quit = true;
            }
            else if (stat == VALIDSETCMD || stat == VALIDZONECMD)
            {
                bool ok = doCommand
                (
                    mesh,
                    setType,
                    setName,
                    actionName,
                    writeVTK,
                    loop,   // if in looping mode dump sets to time directory
                    noSync,
                    is
                );

                if (!ok && batch)
                {
                    // Exit with error.
                    quit = true;
                    status = 1;
                }
            }

        } while (!quit);

        if (quit)
        {
            break;
        }
    }

    Info<< "End\n" << endl;

    return status;
}


// ************************************************************************* //
