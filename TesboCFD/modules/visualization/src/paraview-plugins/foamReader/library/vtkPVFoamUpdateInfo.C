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

\*---------------------------------------------------------------------------*/

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "Cloud.H"
#include "areaFaMesh.H"
#include "OSspecific.H"  // For isFile

// VTK includes
#include "vtkDataArraySelection.h"

// OpenFOAM/VTK interface
#include "vtkPVFoamReader.h"

// Templates (only needed here)
#include "vtkPVFoamUpdateTemplates.C"

// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{
public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(word("regIOobject")))
        {
            close();
        }

    // Member Functions

        bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }
};

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType>
Foam::wordList Foam::vtkPVFoam::getZoneNames
(
    const ZoneMesh<ZoneType, polyMesh>& zmesh
)
{
    wordList names(zmesh.size());
    label nZone = 0;

    for (const auto& zn : zmesh)
    {
        if (!zn.empty())
        {
            names[nZone++] = zn.name();
        }
    }
    names.resize(nZone);

    return names;
}


Foam::wordList Foam::vtkPVFoam::getZoneNames(const word& zoneType) const
{
    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            meshDir_,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    wordList names;
    if (ioObj.typeHeaderOk<cellZoneMesh>(false, false))
    {
        zonesEntries zones(ioObj);

        names.resize(zones.size());
        label nZone = 0;

        for (const auto& zn : zones)
        {
            names[nZone++] = zn.keyword();
        }
    }

    return names;
}


void Foam::vtkPVFoam::updateInfoInternalMesh
(
    vtkDataArraySelection* select
)
{
    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    // Determine mesh parts (internalMesh, patches...)
    // Add internal mesh as first entry
    rangeVolume_.reset(select->GetNumberOfArrays(), 1);
    select->AddArray("internalMesh");

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoAreaMesh
(
    vtkDataArraySelection* select
)
{
    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    rangeArea_.reset(select->GetNumberOfArrays(), 0);

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName faMeshPrefix(faMesh::meshSubDir);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        faMeshPrefix = meshRegion_/faMeshPrefix;
    }

    // Check for finiteArea mesh
    if
    (
        Foam::isFile
        (
            dbPtr_->constantPath()/faMeshPrefix/"faceLabels"
        )
    )
    {
        rangeArea_ += 1;
        select->AddArray("areaMesh");
    }

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoLagrangian
(
    vtkDataArraySelection* select
)
{
    DebugInfo
        << "<beg> " << FUNCTION_NAME << nl
        << "    " << dbPtr_->timePath()/cloud::prefix << nl;

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    // List of lagrangian objects across all times
    HashSet<fileName> names;

    // Get times list. Flush first to force refresh.
    fileHandler().flush();

    for (const instant& t : dbPtr_().times())
    {
        names.insert
        (
            fileHandler().readDir
            (
                dbPtr_->path()/t.name()/lagrangianPrefix,
                fileName::DIRECTORY
            )
        );
    }

    rangeClouds_.reset(select->GetNumberOfArrays());
    rangeClouds_ += addToArray(select, "lagrangian/", names.sortedToc());

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoPatches
(
    vtkDataArraySelection* select,
    HashSet<string>& enabledEntries
)
{
    DebugInfo
        << "<beg> " << FUNCTION_NAME
        << " [volMeshPtr=" << (volMeshPtr_ ? "set" : "null") << "]" << nl;

    rangePatches_.reset(select->GetNumberOfArrays());

    if (volMeshPtr_)
    {
        const polyBoundaryMesh& patches = volMeshPtr_->boundaryMesh();
        const HashTable<labelList>& groups = patches.groupPatchIDs();
        DynamicList<string> displayNames(groups.size());

        // Add (non-zero) patch groups to the list of mesh parts
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAllConstIters(groups, iter)
        {
            const auto& groupName = iter.key();
            const auto& patchIDs  = iter.val();

            label nFaces = 0;
            for (auto patchId : patchIDs)
            {
                nFaces += patches[patchId].size();
            }

            if (!nFaces)
            {
                // Skip if group has no faces
                continue;
            }

            // Valid patch if nFace > 0 - add patch to GUI list
            const string dpyName = "group/" + groupName;
            displayNames.append(dpyName);

            // Optionally replace group with patch name selections
            // - must remove the group from the select itself, otherwise
            //   it can toggle on, but not toggle off very well
            if
            (
                !reader_->GetShowGroupsOnly()
             && enabledEntries.erase(dpyName)
            )
            {
                for (auto patchId : patchIDs)
                {
                    const polyPatch& pp = patches[patchId];
                    if (pp.size())
                    {
                        enabledEntries.insert("patch/" + pp.name());
                    }
                }
            }
        }

        Foam::sort(displayNames);  // Sorted group names
        rangePatches_ += addToArray(select, displayNames);

        // Add (non-zero) patches to the list of mesh parts
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (!reader_->GetShowGroupsOnly())
        {
            for (const polyPatch& pp : patches)
            {
                if (pp.size())
                {
                    // Add patch to GUI list
                    select->AddArray(("patch/" + pp.name()).c_str());
                    ++rangePatches_;
                }
            }
        }
    }
    else
    {
        // mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        IOobject ioObj
        (
            "boundary",
            dbPtr_().findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            dbPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // This should only ever fail if the mesh region doesn't exist
        if (ioObj.typeHeaderOk<polyBoundaryMesh>(true, false))
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            // Read patches, determine sizes and patch groups
            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());
            HashTable<labelHashSet> groups(2*patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();
                wordList groupNames;

                sizes[patchi] = patchDict.get<label>("nFaces");
                names[patchi] = patchEntries[patchi].keyword();

                if
                (
                    sizes[patchi]  // Valid patch if nFace > 0
                 && patchDict.readIfPresent("inGroups", groupNames)
                )
                {
                    for (const auto& groupName : groupNames)
                    {
                        groups(groupName).insert(patchi);
                    }
                }
            }


            // Add (non-zero) patch groups to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            DynamicList<string> displayNames(groups.size());

            forAllConstIters(groups, iter)
            {
                const auto& groupName = iter.key();
                const auto& patchIDs  = iter.val();

                const string dpyName = "group/" + groupName;
                displayNames.append(dpyName);

                // Optionally replace group with patch name selections
                // - must remove the group from the select itself, otherwise
                //   it can toggle on, but not toggle off very well
                if
                (
                    !reader_->GetShowGroupsOnly()
                 && enabledEntries.erase(dpyName)
                )
                {
                    for (auto patchId : patchIDs)
                    {
                        enabledEntries.insert("patch/" + names[patchId]);
                    }
                }
            }

            Foam::sort(displayNames);  // Sorted group names
            rangePatches_ += addToArray(select, displayNames);

            // Add (non-zero) patches to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (!reader_->GetShowGroupsOnly())
            {
                forAll(names, patchi)
                {
                    // Valid patch if nFace > 0 - add patch to GUI list
                    if (sizes[patchi])
                    {
                        select->AddArray(("patch/" + names[patchi]).c_str());
                        ++rangePatches_;
                    }
                }
            }
        }
    }

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoZones
(
    vtkDataArraySelection* select
)
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    DebugInfo
        << "<beg> " << FUNCTION_NAME
        << " [volMeshPtr=" << (volMeshPtr_ ? "set" : "null") << "]" << nl;

    // cellZones
    {
        const wordList names =
        (
            volMeshPtr_
          ? getZoneNames(volMeshPtr_->cellZones())
          : getZoneNames("cellZones")
        );

        rangeCellZones_.reset(select->GetNumberOfArrays());
        rangeCellZones_ += addToArray(select, "cellZone/", names);
    }

    // faceZones
    {
        const wordList names =
        (
            volMeshPtr_
          ? getZoneNames(volMeshPtr_->faceZones())
          : getZoneNames("faceZones")
        );

        rangeFaceZones_.reset(select->GetNumberOfArrays());
        rangeFaceZones_ += addToArray(select, "faceZone/", names);
    }

    // pointZones
    {
        const wordList names =
        (
            volMeshPtr_
          ? getZoneNames(volMeshPtr_->pointZones())
          : getZoneNames("pointZones")
        );

        rangePointZones_.reset(select->GetNumberOfArrays());
        rangePointZones_ += addToArray(select, "pointZone/", names);
    }

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoSets
(
    vtkDataArraySelection* select
)
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    // Add names of sets. Search for last time directory with a sets
    // subdirectory. Take care not to search beyond the last mesh.

    const word facesInstance = dbPtr_().findInstance
    (
        meshDir_,
        "faces",
        IOobject::READ_IF_PRESENT
    );

    const word setsInstance = dbPtr_().findInstance
    (
        meshDir_/"sets",
        word::null,
        IOobject::READ_IF_PRESENT,
        facesInstance
    );

    const IOobjectList objects(dbPtr_(), setsInstance, meshDir_/"sets");

    DebugInfo
        << "     updateInfoSets read "
        << objects.names() << " from " << setsInstance << nl;


    rangeCellSets_.reset(select->GetNumberOfArrays());
    rangeCellSets_ += addToSelection<cellSet>
    (
        select,
        "cellSet/",
        objects
    );

    rangeFaceSets_.reset(select->GetNumberOfArrays());
    rangeFaceSets_ += addToSelection<faceSet>
    (
        select,
        "faceSet/",
        objects
    );

    rangePointSets_.reset(select->GetNumberOfArrays());
    rangePointSets_ += addToSelection<pointSet>
    (
        select,
        "pointSet/",
        objects
    );

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoContinuumFields
(
    vtkDataArraySelection* select
)
{
    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    // First call: enable all
    const bool enableAll = (!select->GetNumberOfArrays() && !volMeshPtr_);

    // Preserve the enabled selections
    HashSet<string> enabled;
    if (!enableAll)
    {
        enabled = getSelectedArraySet(select);
    }
    select->RemoveAllArrays();   // Clear existing list

    // Use the db directly since this might be called without a mesh,
    // but the region name is still required
    word regionPrefix;
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        regionPrefix = meshRegion_;
    }

    // Objects for this time and mesh region
    IOobjectList objects(dbPtr_(), dbPtr_().timeName(), regionPrefix);

    updateInfoFields<fvPatchField, volMesh>(select, objects);
    updateInfoFields<faPatchField, areaMesh>(select, objects);

    if (enableAll)
    {
        select->EnableAllArrays();
    }
    else
    {
        setSelectedArrayEntries(select, enabled);
    }

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoPointFields
(
    vtkDataArraySelection* select
)
{
    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    // Preserve the enabled selections
    HashSet<string> enabled = getSelectedArraySet(select);

    select->RemoveAllArrays();   // Clear existing list

    // Use the db directly since this might be called without a mesh,
    // but the region name is still required
    word regionPrefix;
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        regionPrefix = meshRegion_;
    }

    // Objects for this time and mesh region
    IOobjectList objects(dbPtr_(), dbPtr_().timeName(), regionPrefix);

    updateInfoFields<pointPatchField, pointMesh>(select, objects);

    setSelectedArrayEntries(select, enabled);  // Adjust/restore selected

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


void Foam::vtkPVFoam::updateInfoLagrangianFields
(
    vtkDataArraySelection* select
)
{
    DebugInfo << "<beg> " << FUNCTION_NAME << nl;

    // Preserve the enabled selections
    HashSet<string> enabled = getSelectedArraySet(select);
    select->RemoveAllArrays();

    const arrayRange& range = rangeClouds_;
    if (range.empty())
    {
        return;
    }

    // Reuse the previously determined cloud information.
    DynamicList<word> cloudNames(range.size());
    for (auto partId : range)
    {
        cloudNames.append(getReaderPartName(partId));
    }

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    // List of lagrangian fields across all clouds and all times.
    // ParaView displays "(partial)" after field names that only apply
    // to some of the clouds.
    HashTable<wordHashSet> fields;

    fileHandler().flush();
    for (const instant& t : dbPtr_().times())
    {
        for (const auto& cloudName : cloudNames)
        {
            const HashTable<wordHashSet> localFields =
                IOobjectList
                (
                    dbPtr_(),
                    t.name(),
                    lagrangianPrefix/cloudName
                ).classes();

            forAllConstIters(localFields, iter)
            {
                fields(iter.key()) |= iter.val();
            }
        }
    }

    // Known/supported field-types
    addToSelection<IOField<label>>(select, fields);
    addToSelection<IOField<scalar>>(select, fields);
    addToSelection<IOField<vector>>(select, fields);
    addToSelection<IOField<sphericalTensor>>(select, fields);
    addToSelection<IOField<symmTensor>>(select, fields);
    addToSelection<IOField<tensor>>(select, fields);

    // Restore the enabled selections
    setSelectedArrayEntries(select, enabled);

    DebugInfo << "<end> " << FUNCTION_NAME << nl;
}


// ************************************************************************* //
