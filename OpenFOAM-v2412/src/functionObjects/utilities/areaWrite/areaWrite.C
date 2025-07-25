/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "areaWrite.H"
#include "polySurface.H"

#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "areaFields.H"
#include "HashOps.H"
#include "ListOps.H"
#include "Time.H"
#include "IndirectList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(areaWrite, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        areaWrite,
        dictionary
    );
}

Foam::scalar Foam::areaWrite::mergeTol_ = 1e-10;


// Implementation
#include "areaWriteImpl.C"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::areaWrite::areaWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    loadFromFiles_(false),
    verbose_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    selectAreas_(),
    fieldSelection_(),
    meshes_(),
    surfaces_(nullptr),
    writers_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


Foam::areaWrite::areaWrite
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::fvMeshFunctionObject(name, obr, dict),
    loadFromFiles_(loadFromFiles),
    verbose_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    selectAreas_(),
    fieldSelection_(),
    meshes_(),
    surfaces_(nullptr),
    writers_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::areaWrite::verbose(const bool on) noexcept
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


bool Foam::areaWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    writers_.clear();
    selectAreas_.clear();
    fieldSelection_.clear();

    surfaces_.reset
    (
        new objectRegistry
        (
            IOobject
            (
                "::areaWrite::",
                obr_.time().constant(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );

    verbose_ = dict.getOrDefault("verbose", false);

    // Registry containing all finite-area meshes on the polyMesh
    const auto* faRegistry = faMesh::registry(mesh_);

    dict.readIfPresent("areas", selectAreas_);

    if (selectAreas_.empty())
    {
        word areaName;

        if (!dict.readIfPresent("area", areaName))
        {
            if (faRegistry)
            {
                wordList available = faRegistry->sortedNames<faMesh>();
                if (!available.empty())
                {
                    areaName = available.front();
                }
            }
        }

        if (!areaName.empty())
        {
            selectAreas_.resize(1);
            selectAreas_.front() = areaName;
        }
    }

    // Restrict to specified meshes
    meshes_.clear();

    if (faRegistry)
    {
        meshes_ = faRegistry->csorted<faMesh>(selectAreas_);
    }

    dict.readEntry("fields", fieldSelection_);
    fieldSelection_.uniq();


    // Surface writer type and format options
    const word writerType = dict.get<word>("surfaceFormat");

    const dictionary writerOptions
    (
        surfaceWriter::formatOptions(dict, writerType)
    );

    for (const faMesh& areaMesh : meshes_)
    {
        const word& areaName = areaMesh.name();

        // Define surface writer, but do NOT yet attach a surface

        auto surfWriter = surfaceWriter::New(writerType, writerOptions);

        // Use outputDir/TIME/surface-name
        surfWriter->useTimeDir(true);
        surfWriter->verbose(verbose_);

        writers_.set(areaName, surfWriter);
    }

    // Ensure all surfaces and merge information are expired
    expire();

    return true;
}


bool Foam::areaWrite::execute()
{
    return true;
}


bool Foam::areaWrite::write()
{
    // Just needed for warnings
    wordList allFields;
    HashTable<wordHashSet> selected;
    DynamicList<label> missed(fieldSelection_.size());


    for (const faMesh& areaMesh : meshes_)
    {
        const word& areaName = areaMesh.name();

        polySurface* surfptr = surfaces_->getObjectPtr<polySurface>(areaName);

        if (!surfptr)
        {
            // Construct null and add to registry (owned by registry)
            surfptr = new polySurface(areaName, *surfaces_, true);
        }

        pointField pts(areaMesh.patch().localPoints());
        faceList fcs(areaMesh.patch().localFaces());

        // Copy in geometry
        surfptr->transfer(std::move(pts), std::move(fcs));

        surfaceWriter& outWriter = *writers_[areaName];

        if (outWriter.needsUpdate())
        {
            outWriter.setSurface(*surfptr);
        }


        // Determine the per-surface number of fields
        // Only seems to be needed for VTK legacy

        selected.clear();

        IOobjectList objects;

        if (loadFromFiles_)
        {
            // Check files for a particular time
            objects = IOobjectList(areaMesh.thisDb(), obr_.time().timeName());

            allFields = objects.names();
            selected = objects.classes(fieldSelection_);
        }
        else
        {
            // Check currently available fields
            allFields = areaMesh.thisDb().names();
            selected = areaMesh.thisDb().classes(fieldSelection_);
        }

        // Parallel consistency (no-op in serial)
        Pstream::mapCombineReduce(selected, HashSetOps::plusEqOp<word>());

        missed.clear();

        // Detect missing fields
        forAll(fieldSelection_, i)
        {
            if (!ListOps::found(allFields, fieldSelection_[i]))
            {
                missed.push_back(i);
            }
        }

        if (missed.size())
        {
            WarningInFunction
                << nl
                << "Cannot find "
                << (loadFromFiles_ ? "field file" : "registered field")
                << " matching "
                << UIndirectList<wordRe>(fieldSelection_, missed) << endl;
        }


        // Currently only support area field types
        label nAreaFields = 0;

        forAllConstIters(selected, iter)
        {
            const word& clsName = iter.key();
            const label n = iter.val().size();

            if
            (
                fieldTypes::area.contains(clsName)
             || fieldTypes::area_internal.contains(clsName)
            )
            {
                nAreaFields += n;
            }
        }


        // Propagate field counts (per surface)
        outWriter.nFields(nAreaFields);


        // Begin writing

        outWriter.open(outputPath_/areaName);

        outWriter.beginTime(obr_.time());

        // Write fields

        {
            // Area fields
            #undef  doLocalCode
            #define doLocalCode(Type)                                         \
            performAction                                                     \
            <                                                                 \
                GeometricField<Type, Foam::faPatchField, Foam::areaMesh>      \
            >                                                                 \
            (                                                                 \
                outWriter, areaMesh, objects                                  \
            );                                                                \

            doLocalCode(scalar);
            doLocalCode(vector);
            doLocalCode(sphericalTensor);
            doLocalCode(symmTensor);
            doLocalCode(tensor);

            // Area internal fields
            #undef  doLocalCode
            #define doLocalCode(Type)                                         \
            performAction                                                     \
            <                                                                 \
                DimensionedField<Type, Foam::areaMesh>                        \
            >                                                                 \
            (                                                                 \
                outWriter, areaMesh, objects                                  \
            );

            doLocalCode(scalar);
            doLocalCode(vector);
            doLocalCode(sphericalTensor);
            doLocalCode(symmTensor);
            doLocalCode(tensor);

            #undef doLocalCode
        }

        // Finish this time step

        // Write geometry if no fields were written so that we still
        // can have something to look at

        if (!outWriter.wroteData())
        {
            outWriter.write();
        }

        outWriter.endTime();
    }

    return true;
}


void Foam::areaWrite::expire()
{
    surfaces_->clear();

    // Dimension as fraction of mesh bounding box
    const scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    forAllIters(writers_, iter)
    {
        surfaceWriter& writer = *(iter.val());
        writer.expire();
        writer.mergeDim(mergeDim);
    }
}


void Foam::areaWrite::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }
}


void Foam::areaWrite::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::areaWrite::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


Foam::scalar Foam::areaWrite::mergeTol() noexcept
{
    return mergeTol_;
}


Foam::scalar Foam::areaWrite::mergeTol(const scalar tol) noexcept
{
    scalar old(mergeTol_);
    mergeTol_ = tol;
    return old;
}


// ************************************************************************* //
