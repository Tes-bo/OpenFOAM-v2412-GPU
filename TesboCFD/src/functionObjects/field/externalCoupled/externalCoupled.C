/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "externalCoupled.H"
#include "stringListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"
#include "Fstream.H"
#include "SpanStream.H"
#include "volFields.H"
#include "globalIndex.H"
#include "fvMesh.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(externalCoupled, 0);
    addToRunTimeSelectionTable(functionObject, externalCoupled, dictionary);
}
}

Foam::string Foam::functionObjects::externalCoupled::patchKey = "// Patch:";


namespace Foam
{
//! \cond fileScope
//- Write list content with size, bracket, content, bracket one-per-line.
//  This makes for consistent for parsing, regardless of the list length.
template <class T>
static void writeList(Ostream& os, const string& header, const UList<T>& L)
{
    // Header string
    os  << header.c_str() << nl;

    // Write size and start delimiter
    os  << L.size() << nl
        << token::BEGIN_LIST << nl;

    // Write contents
    forAll(L, i)
    {
        os << L[i] << nl;
    }

    // Write end delimiter
    os << token::END_LIST << nl << endl;
}
//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::externalCoupled::groupDir
(
    const fileName& commsDir,
    const word& regionGroupName,
    const wordRe& groupName
)
{
    fileName result
    (
        commsDir
      / regionGroupName
      / word::validate(groupName)
    );
    result.clean();  // Remove unneeded ".."

    return result;
}


void Foam::functionObjects::externalCoupled::readColumns
(
    const label nRows,
    const label nColumns,
    autoPtr<IFstream>& masterFilePtr,
    List<scalarField>& data
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(globalIndex::gatherOnly{}, nRows);

    PstreamBuffers pBufs;

    if (UPstream::master())
    {
        // Read data from file and send to destination processor
        auto& ifile = masterFilePtr();
        string line;

        for (const int proci : UPstream::allProcs())
        {
            // Temporary storage
            List<scalarField> values(nColumns);

            // Number of rows to read for processor proci
            const label procNRows = globalFaces.localSize(proci);

            forAll(values, columni)
            {
                values[columni].setSize(procNRows);
            }

            for (label rowi = 0; rowi < procNRows; ++rowi)
            {
                // Get a line
                do
                {
                    if (!ifile.good())
                    {
                        FatalIOErrorInFunction(ifile)
                            << "Trying to read data for processor " << proci
                            << " row " << rowi
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.totalSize()
                            << ") ?" << exit(FatalIOError);
                    }

                    ifile.getLine(line);
                }
                while (line.empty() || line[0] == '#');

                ISpanStream isstr(line);

                for (label columni = 0; columni < nColumns; ++columni)
                {
                    isstr >> values[columni][rowi];
                }
            }

            // Send to proci
            UOPstream toProc(proci, pBufs);
            toProc << values;
        }
    }
    pBufs.finishedScatters();

    // Get scattered data
    UIPstream fromMaster(UPstream::masterNo(), pBufs);
    fromMaster >> data;
}


void Foam::functionObjects::externalCoupled::readLines
(
    const label nRows,
    autoPtr<IFstream>& masterFilePtr,
    std::string& lines
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(globalIndex::gatherOnly{}, nRows);

    PstreamBuffers pBufs;

    if (UPstream::master())
    {
        // Read line from file and send to destination processor
        auto& ifile = masterFilePtr();
        string line;

        for (const int proci : UPstream::allProcs())
        {
            // Number of rows to read for processor proci
            const label procNRows = globalFaces.localSize(proci);

            UOPstream toProc(proci, pBufs);

            for (label rowi = 0; rowi < procNRows; ++rowi)
            {
                // Get a line
                do
                {
                    if (!ifile.good())
                    {
                        FatalIOErrorInFunction(ifile)
                            << "Trying to read data for processor " << proci
                            << " row " << rowi
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.totalSize()
                            << ") ?" << exit(FatalIOError);
                    }

                    ifile.getLine(line);
                }
                while (line.empty() || line[0] == '#');

                // Send line (without newline) to the destination processor
                toProc << line;
            }
        }
    }

    pBufs.finishedScatters();

    // Sizing is approximate (slightly too large)
    lines.reserve(pBufs.recvDataCount(UPstream::masterNo()));

    // Get scattered data
    UIPstream fromMaster(UPstream::masterNo(), pBufs);
    for (label rowi = 0; rowi < nRows; ++rowi)
    {
        string line(fromMaster);
        lines += line;
        lines += '\n';
    }
}


void Foam::functionObjects::externalCoupled::writeGeometry
(
    const UPtrList<const fvMesh>& meshes,
    const fileName& commsDir,
    const wordRe& groupName
)
{
    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // Make sure meshes are provided in sorted order
    checkOrder(regionNames);

    fileName dir(groupDir(commsDir, compositeName(regionNames), groupName));

    autoPtr<OFstream> osPointsPtr;
    autoPtr<OFstream> osFacesPtr;
    if (UPstream::master())
    {
        Foam::mkDir(dir);
        osPointsPtr.reset(new OFstream(dir/"patchPoints"));
        osFacesPtr.reset(new OFstream(dir/"patchFaces"));

        osPointsPtr() << "// Group: " << groupName << endl;
        osFacesPtr()  << "// Group: " << groupName << endl;

        Info<< typeName << ": writing geometry to " << dir << endl;
    }

    // Individual region/patch entries

    DynamicList<face> allFaces;
    DynamicField<point> allPoints;

    labelList pointToGlobal;
    labelList uniquePointIDs;
    for (const fvMesh& mesh : meshes)
    {
        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet
            (
                wordRes(Foam::one{}, groupName)
            ).sortedToc()
        );

        for (const label patchi : patchIDs)
        {
            const polyPatch& p = mesh.boundaryMesh()[patchi];

            mesh.globalData().mergePoints
            (
                p.meshPoints(),
                p.meshPointMap(),
                pointToGlobal,
                uniquePointIDs
            );

            label proci = Pstream::myProcNo();

            List<pointField> collectedPoints(Pstream::nProcs());
            collectedPoints[proci] = pointField(mesh.points(), uniquePointIDs);
            Pstream::gatherList(collectedPoints);

            List<faceList> collectedFaces(Pstream::nProcs());
            faceList& patchFaces = collectedFaces[proci];
            patchFaces = p.localFaces();
            for (auto& f : patchFaces)
            {
                inplaceRenumber(pointToGlobal, f);
            }
            Pstream::gatherList(collectedFaces);

            if (UPstream::master())
            {
                allPoints.clear();
                allFaces.clear();

                for (const int proci : Pstream::allProcs())
                {
                    allPoints.append(collectedPoints[proci]);
                    allFaces.append(collectedFaces[proci]);
                }

                Info<< typeName << ": mesh " << mesh.name()
                    << ", patch " << p.name()
                    << ": writing " << allPoints.size() << " points to "
                    << osPointsPtr().name() << nl
                    << typeName << ": mesh " << mesh.name()
                    << ", patch " << p.name()
                    << ": writing " << allFaces.size() << " faces to "
                    << osFacesPtr().name() << endl;

                // The entry name (region / patch)
                const string entryHeader =
                    patchKey + ' ' + mesh.name() + ' ' + p.name();

                writeList(osPointsPtr(), entryHeader, allPoints);
                writeList(osFacesPtr(),  entryHeader, allFaces);
            }
        }
    }
}


Foam::word Foam::functionObjects::externalCoupled::compositeName
(
    const wordList& regionNames
)
{
    if (regionNames.size() == 0)
    {
        FatalErrorInFunction
            << "Empty regionNames" << abort(FatalError);
        return word::null;
    }
    else if (regionNames.size() == 1)
    {
        // For compatibility with single region cases
        // - suppress single region name
        return polyMesh::regionName(regionNames[0]);
    }

    // Enforce lexical ordering
    checkOrder(regionNames);

    word composite(regionNames[0]);
    for (label i = 1; i < regionNames.size(); ++i)
    {
        composite += "_" + regionNames[i];
    }

    return composite;
}


void Foam::functionObjects::externalCoupled::checkOrder
(
    const wordList& regionNames
)
{
    labelList order(sortedOrder(regionNames));
    if (order != identity(regionNames.size()))
    {
        FatalErrorInFunction
            << "regionNames " << regionNames << " not in alphabetical order :"
            << order << exit(FatalError);
    }
}


void Foam::functionObjects::externalCoupled::initCoupling()
{
    if (initialisedCoupling_)
    {
        return;
    }

    // Write the geometry if not already there
    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];
        const wordList& regionNames = regionGroupRegions_[regioni];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, regi)
        {
            meshes.set(regi, time_.findObject<fvMesh>(regionNames[regi]));
        }

        const labelList& groups = regionToGroups_[compName];

        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];

            bool geomExists = false;
            if (Pstream::master())
            {
                fileName dir(groupDir(commDirectory(), compName, groupName));

                geomExists =
                    isFile(dir/"patchPoints")
                 || isFile(dir/"patchFaces");
            }

            Pstream::broadcast(geomExists);

            if (!geomExists)
            {
                writeGeometry(meshes, commDirectory(), groupName);
            }
        }
    }

    if (slaveFirst())
    {
        // Wait for initial data to be made available
        waitForSlave();

        // Read data passed back from external source
        readDataMaster();
    }

    initialisedCoupling_ = true;
}


void Foam::functionObjects::externalCoupled::performCoupling()
{
    // Ensure coupling has been initialised
    initCoupling();

    // Write data for external source
    writeDataMaster();

    // Signal external source to execute (by removing lock file)
    // - Wait for slave to provide data
    useSlave();

    // Wait for response - and catch any abort information sent from slave
    const auto action = waitForSlave();

    // Remove old data files from OpenFOAM
    removeDataMaster();

    // Read data passed back from external source
    readDataMaster();

    // Signal external source to wait (by creating the lock file)
    useMaster();

    // Update information about last triggering
    lastTrigger_ = time_.timeIndex();

    // Process any abort information sent from slave
    if
    (
        action != time_.stopAt()
     && action != Time::stopAtControls::saUnknown
    )
    {
        Info<< type() << ": slave requested action "
            << Time::stopAtControlNames[action] << endl;

        time_.stopAt(action);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::externalCoupled::externalCoupled
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    externalFileCoupler(),
    calcFrequency_(-1),
    lastTrigger_(-1),
    initialisedCoupling_(false)
{
    read(dict);

    if (!slaveFirst())
    {
        useMaster();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::externalCoupled::execute()
{
    // Not initialized or overdue
    if
    (
        !initialisedCoupling_
     || (time_.timeIndex() >= lastTrigger_ + calcFrequency_)
    )
    {
        performCoupling();
    }

    return false;
}


bool Foam::functionObjects::externalCoupled::execute(const label subIndex)
{
    performCoupling();

    return true;
}


bool Foam::functionObjects::externalCoupled::end()
{
    functionObject::end();

    // Remove old data files
    removeDataMaster();
    removeDataSlave();
    shutdown();

    return true;
}


bool Foam::functionObjects::externalCoupled::read(const dictionary& dict)
{
    timeFunctionObject::read(dict);
    externalFileCoupler::readDict(dict);

    calcFrequency_ =
        dict.getCheckOrDefault("calcFrequency", 1, labelMinMax::ge(1));

    // Leave trigger intact

    // Get names of all fvMeshes (and derived types)
    wordList allRegionNames(time_.sortedNames<fvMesh>());

    const dictionary& allRegionsDict = dict.subDict("regions");
    for (const entry& dEntry : allRegionsDict)
    {
        if (!dEntry.isDict())
        {
            FatalIOErrorInFunction(allRegionsDict)
                << "Regions must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const wordRe regionGroupName(dEntry.keyword());
        const dictionary& regionDict = dEntry.dict();

        labelList regionIDs
        (
            wordRes::matching(regionGroupName, allRegionNames)
        );

        const wordList regionNames(allRegionNames, regionIDs);

        regionGroupNames_.append(compositeName(regionNames));
        regionGroupRegions_.append(regionNames);

        for (const entry& dEntry : regionDict)
        {
            if (!dEntry.isDict())
            {
                FatalIOErrorInFunction(regionDict)
                    << "Regions must be specified in dictionary format"
                    << exit(FatalIOError);
            }

            const wordRe groupName(dEntry.keyword());
            const dictionary& groupDict = dEntry.dict();

            const label nGroups = groupNames_.size();
            const wordList readFields(groupDict.get<wordList>("readFields"));
            const wordList writeFields(groupDict.get<wordList>("writeFields"));

            auto fnd = regionToGroups_.find(regionGroupNames_.last());
            if (fnd.good())
            {
                fnd().append(nGroups);
            }
            else
            {
                regionToGroups_.insert
                (
                    regionGroupNames_.last(),
                    labelList(one{}, nGroups)
                );
            }
            groupNames_.append(groupName);
            groupReadFields_.append(readFields);
            groupWriteFields_.append(writeFields);
        }
    }


    Info<< type() << ": Communicating with regions:" << endl;
    for (const word& compName : regionGroupNames_)
    {
        Info<< "Region: " << compName << nl << incrIndent;
        const labelList& groups = regionToGroups_[compName];
        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];

            Info<< indent << "patchGroup: " << groupName << "\t"
                << nl
                << incrIndent
                << indent << "Reading fields: "
                << groupReadFields_[groupi]
                << nl
                << indent << "Writing fields: "
                << groupWriteFields_[groupi]
                << nl
                << decrIndent;
        }
        Info<< decrIndent;
    }
    Info<< endl;


    // Note: we should not have to make directories since the geometry
    //       should already be written - but just make sure
    if (Pstream::master())
    {
        for (const word& compName : regionGroupNames_)
        {
            const labelList& groups = regionToGroups_[compName];
            for (const label groupi : groups)
            {
                const wordRe& groupName = groupNames_[groupi];

                fileName dir(groupDir(commDirectory(), compName, groupName));

                if (!isDir(dir))
                {
                    Log << type() << ": creating communications directory "
                        << dir << endl;
                    mkDir(dir);
                }
            }
        }
    }

    return true;
}


void Foam::functionObjects::externalCoupled::readDataMaster()
{
    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];
        const wordList& regionNames = regionGroupRegions_[regioni];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, regi)
        {
            meshes.set(regi, time_.findObject<fvMesh>(regionNames[regi]));
        }

        const labelList& groups = regionToGroups_[compName];

        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupReadFields_[groupi];

            for (const word& fieldName : fieldNames)
            {
                const bool ok =
                (
                    readData<scalar>(meshes, groupName, fieldName)
                 || readData<vector>(meshes, groupName, fieldName)
                 || readData<sphericalTensor>(meshes, groupName, fieldName)
                 || readData<symmTensor>(meshes, groupName, fieldName)
                 || readData<tensor>(meshes, groupName, fieldName)
                );

                if (!ok)
                {
                    WarningInFunction
                        << "Field " << fieldName << " in regions " << compName
                        << " was not found." << endl;
                }
            }
        }
    }
}


void Foam::functionObjects::externalCoupled::writeDataMaster() const
{
    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];
        const wordList& regionNames = regionGroupRegions_[regioni];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, regi)
        {
            meshes.set(regi, time_.findObject<fvMesh>(regionNames[regi]));
        }

        const labelList& groups = regionToGroups_[compName];

        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupWriteFields_[groupi];

            for (const word& fieldName : fieldNames)
            {
                const bool ok =
                (
                    writeData<scalar>(meshes, groupName, fieldName)
                 || writeData<vector>(meshes, groupName, fieldName)
                 || writeData<sphericalTensor>(meshes, groupName, fieldName)
                 || writeData<symmTensor>(meshes, groupName, fieldName)
                 || writeData<tensor>(meshes, groupName, fieldName)
                );

                if (!ok)
                {
                    WarningInFunction
                        << "Field " << fieldName << " in regions " << compName
                        << " was not found." << endl;
                }
            }
        }
    }
}


void Foam::functionObjects::externalCoupled::removeDataMaster() const
{
    if (!Pstream::master())
    {
        return;
    }

    Log << type() << ": removing data files written by master" << nl;

    for (const word& compName : regionGroupNames_)
    {
        const labelList& groups = regionToGroups_[compName];
        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupReadFields_[groupi];

            for (const word& fieldName : fieldNames)
            {
                Foam::rm
                (
                    groupDir(commDirectory(), compName, groupName)
                  / fieldName + ".out"
                );
            }
        }
    }
}


void Foam::functionObjects::externalCoupled::removeDataSlave() const
{
    if (!Pstream::master())
    {
        return;
    }

    Log << type() << ": removing data files written by slave" << nl;

    for (const word& compName : regionGroupNames_)
    {
        const labelList& groups = regionToGroups_[compName];
        for (const label groupi : groups)
        {
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupReadFields_[groupi];

            for (const word& fieldName : fieldNames)
            {
                Foam::rm
                (
                    groupDir(commDirectory(), compName, groupName)
                  / fieldName + ".in"
                );
            }
        }
    }
}


bool Foam::functionObjects::externalCoupled::write()
{
    return true;
}


// ************************************************************************* //
