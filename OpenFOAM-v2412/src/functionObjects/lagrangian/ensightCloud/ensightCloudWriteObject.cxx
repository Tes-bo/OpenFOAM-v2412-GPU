/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "ensightCloudWriteObject.H"
#include "ensightCells.H"
#include "Cloud.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "ensightOutputCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "pointList.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ensightCloudWriteObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ensightCloudWriteObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "ensightCloudWriteObjectImpl.cxx"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ensightCloudWriteObject::writeCloud
(
    const word& cloudName
)
{
    applyFilter_ = false;
    procAddr_.clear();

    const auto* cloudPtr = mesh_.cfindObject<cloud>(cloudName);
    if (!cloudPtr)
    {
        return false;
    }

    const auto& currCloud = *cloudPtr;

    objectRegistry obrTmp
    (
        IOobject
        (
            "ensight::ensightCloud::" + cloudName,
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    currCloud.writeObjects(obrTmp);

    const auto* pointsPtr = cloud::findIOPosition(obrTmp);

    if (!pointsPtr)
    {
        // This should be impossible
        return false;
    }

    applyFilter_ = calculateFilter(obrTmp, log);
    Pstream::reduceOr(applyFilter_);

    // Number of parcels (locally)
    const label nParcels =
    (
        applyFilter_ ? parcelAddr_.count() : pointsPtr->size()
    );

    // Gather sizes (offsets irrelevant)
    procAddr_.reset(globalIndex::gatherOnly{}, nParcels);

    bool noCloud(!procAddr_.totalSize());
    Pstream::broadcast(noCloud);

    if (applyFilter_)
    {
        // Report filtered/unfiltered count
        Log << "After filtering using "
            << procAddr_.totalSize() << '/'
            << (returnReduce(pointsPtr->size(), sumOp<label>()))
            << " parcels" << nl;
    }

    if (pruneEmpty_ && noCloud)
    {
        return false;
    }


    // Copy positions (for simplicity and for filtering).
    // Store as floatVector, since that is what Ensight will write anyhow

    DynamicList<floatVector> positions;
    positions.reserve(UPstream::master() ? procAddr_.maxSize() : nParcels);

    {
        const auto& points = *pointsPtr;

        positions.resize_nocopy(nParcels);

        auto iter = positions.begin();

        if (applyFilter_)
        {
            if (std::is_same<float, vector::cmptType>::value)
            {
                for (const label idx : parcelAddr_)
                {
                    *iter = points[idx];
                    ++iter;
                }
            }
            else
            {
                for (const label idx : parcelAddr_)
                {
                    const auto& pos = points[idx];

                    (*iter).x() = narrowFloat(pos.x());
                    (*iter).y() = narrowFloat(pos.y());
                    (*iter).z() = narrowFloat(pos.z());
                    ++iter;
                }
            }
        }
        else
        {
            if (std::is_same<float, vector::cmptType>::value)
            {
                for (const auto& pos : points)
                {
                    *iter = pos;
                    ++iter;
                }
            }
            else
            {
                for (const auto& pos : points)
                {
                    (*iter).x() = narrowFloat(pos.x());
                    (*iter).y() = narrowFloat(pos.y());
                    (*iter).z() = narrowFloat(pos.z());
                    ++iter;
                }
            }
        }
    }


    // Write positions
    {
        autoPtr<ensightFile> os = ensCase().newCloud(cloudName);

        ensightOutput::writeCloudPositions
        (
            os.ref(),
            positions,
            procAddr_
        );
    }

    // Prevent any possible conversion of positions as a field
    obrTmp.filterKeys
    (
        [](const word& k)
        {
            return k.starts_with("position") || k.starts_with("coordinate");
        },
        true  // prune
    );


    // Write fields

    DynamicList<word> written(obrTmp.size() + currCloud.objectRegistry::size());

    written.push_back
    (
        writeFields<label>(cloudName, obrTmp)
    );
    written.push_back
    (
        writeFields<scalar>(cloudName, obrTmp)
    );
    written.push_back
    (
        writeFields<vector>(cloudName, obrTmp)
    );

    // Any cloudFunctions results
    written.push_back
    (
        writeFields<scalar>(cloudName, currCloud)
    );

    // Record information into the state (all processors)
    //
    // foName
    // {
    //     cloudName
    //     {
    //         file   "<case>/postProcessing/name/casename.case";
    //         fields (U T rho);
    //     }
    // }

    const fileName& file = ensCase().path();

    // Case-local file name with "<case>" to make relocatable
    dictionary propsDict;
    propsDict.add
    (
        "file",
        time_.relativePath(file, true)
    );
    propsDict.add("fields", written);

    setObjectProperty(name(), cloudName, propsDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ensightCloudWriteObject::ensightCloudWriteObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    caseOpts_("format", dict, IOstreamOption::BINARY),
    outputDir_(),
    consecutive_(false),
    pruneEmpty_(false),
    applyFilter_(false),
    procAddr_()
{
    // May still want this?
    // if (postProcess)
    // {
    //     // Disable for post-process mode.
    //     // Emit as FatalError for the try/catch in the caller.
    //     FatalError
    //         << type() << " disabled in post-process mode"
    //         << exit(FatalError);
    // }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightCloudWriteObject::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    // Case/writer options
    consecutive_ = dict.getOrDefault("consecutive", false);

    caseOpts_.width(dict.getOrDefault<label>("width", 8));
    caseOpts_.overwrite(dict.getOrDefault("overwrite", false));

    caseOpts_.timeFormat("timeFormat", dict);
    caseOpts_.timePrecision("timePrecision", dict);


    pruneEmpty_ = dict.getOrDefault("prune", false);

    selectClouds_.clear();
    dict.readIfPresent("clouds", selectClouds_);
    selectClouds_.uniq();
    if (selectClouds_.empty())
    {
        word cloudName;
        if (dict.readIfPresent("cloud", cloudName))
        {
            selectClouds_.push_back(std::move(cloudName));
        }
    }

    selectFields_.clear();
    dict.readIfPresent("fields", selectFields_);
    selectFields_.uniq();

    // Actions to define selection
    parcelSelect_ = dict.subOrEmptyDict("selection");


    // Output directory

    outputDir_.clear();
    dict.readIfPresent("directory", outputDir_);

    if (outputDir_.size())
    {
        // User-defined output directory
        outputDir_.expand();
        if (!outputDir_.isAbsolute())
        {
            outputDir_ = time_.globalPath()/outputDir_;
        }
    }
    else
    {
        // Standard postProcessing/ naming
        outputDir_ = time_.globalPath()/functionObject::outputPrefix/name();
    }
    outputDir_.clean();  // Remove unneeded ".."

    return true;
}


bool Foam::functionObjects::ensightCloudWriteObject::execute()
{
    return true;
}


bool Foam::functionObjects::ensightCloudWriteObject::write()
{
    const wordList cloudNames
    (
        selectClouds_.empty()
      ? mesh_.sortedNames<cloud>()
      : mesh_.sortedNames<cloud>(selectClouds_)
    );

    if (cloudNames.empty())
    {
        return true;  // skip - nothing available
    }

    if (!ensCase_)
    {
        ensCase_.reset
        (
            new ensightCase(outputDir_, time_.globalCaseName(), caseOpts_)
        );

        // Generate a (non-moving) dummy geometry
        // - ParaView ensight-reader needs this, and usually ensight does too
        autoPtr<ensightGeoFile> os = ensCase().newGeometry(false);

        if (os)
        {
            os->beginGeometry();
        }
        ensightCells::writeBox(os.ref(), mesh_.bounds());
    }

    if (consecutive_)
    {
        ensCase().nextTime(time_.value());
    }
    else
    {
        ensCase().setTime(time_.value(), time_.timeIndex());
    }

    Log << type() << ' ' << name() << " write" << nl;

    // Each cloud separately
    for (const word& cloudName : cloudNames)
    {
        // writeCloud() includes mkDir (on master)

        if (writeCloud(cloudName))
        {
            Log << "    cloud  : " << endl;
        }
    }

    ensCase().write();  // Flush case information

    return true;
}


// ************************************************************************* //
