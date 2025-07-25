/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "dataCloud.H"
#include "Cloud.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "pointList.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dataCloud, 0);
    addToRunTimeSelectionTable(functionObject, dataCloud, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::dataCloud::writeCloud
(
    const fileName& outputName,
    const word& cloudName
)
{
    applyFilter_ = false;

    const auto* cloudPtr = mesh_.findObject<cloud>(cloudName);
    if (!cloudPtr)
    {
        return false;
    }

    const auto& currCloud = *cloudPtr;

    objectRegistry obrTmp
    (
        IOobject
        (
            "tmp::dataCloud::" + cloudName,
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
    const label nParcels
    (
        applyFilter_ ? parcelAddr_.count() : pointsPtr->size()
    );

    // Total number of parcels on all processes
    const label nTotParcels = returnReduce(nParcels, sumOp<label>());

    if (applyFilter_)
    {
        // Report filtered/unfiltered count
        Log << "After filtering using " << nTotParcels << '/'
            << (returnReduce(pointsPtr->size(), sumOp<label>()))
            << " parcels" << nl;
    }

    if (!nTotParcels)
    {
        return false;
    }

    if (UPstream::master())
    {
        Foam::mkDir(outputName.path());
    }

    return
    (
        writeField<label>(outputName, obrTmp)
     || writeField<scalar>(outputName, obrTmp)
     || writeField<vector>(outputName, obrTmp)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dataCloud::dataCloud
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    printf_(),
    precision_(IOstream::defaultPrecision()),
    applyFilter_(false),
    selectClouds_(),
    fieldName_(),
    directory_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dataCloud::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    const int padWidth = dict.getOrDefault<int>("width", 8);

    // Appropriate printf format - Enforce min/max sanity limits
    if (padWidth < 1 || padWidth > 31)
    {
        printf_.clear();
    }
    else
    {
        printf_ = "%0" + std::to_string(padWidth) + "d";
    }

    precision_ =
        dict.getOrDefault("precision", IOstream::defaultPrecision());


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

    dict.readEntry("field", fieldName_);

    // Actions to define selection
    parcelSelect_ = dict.subOrEmptyDict("selection");

    // Output directory

    directory_.clear();
    dict.readIfPresent("directory", directory_);

    if (directory_.size())
    {
        // User-defined output directory
        directory_.expand();
        if (!directory_.isAbsolute())
        {
            directory_ = time_.globalPath()/directory_;
        }
    }
    else
    {
        // Standard postProcessing/ naming
        directory_ = time_.globalPath()/functionObject::outputPrefix/name();
    }
    directory_.clean();  // Remove unneeded ".."

    return true;
}


bool Foam::functionObjects::dataCloud::execute()
{
    return true;
}


bool Foam::functionObjects::dataCloud::write()
{
    const wordList cloudNames
    (
        selectClouds_.empty()
      ? mesh_.sortedNames<cloud>()
      : mesh_.sortedNames<cloud>(selectClouds_)
    );

    if (cloudNames.empty())
    {
        return true;  // skip - not available
    }

    const word timeDesc = "_" +
    (
        printf_.empty()
      ? Foam::name(time_.timeIndex())
      : word::printf(printf_, time_.timeIndex())
    );

    Log << name() << " output Time: " << time_.timeName() << nl;

    // Each cloud separately
    for (const word& cloudName : cloudNames)
    {
        // Legacy is not to be supported

        const fileName outputName
        (
            directory_/cloudName + timeDesc + ".dat"
        );

        // writeCloud() includes mkDir (on master)

        if (writeCloud(outputName, cloudName))
        {
            Log << "    cloud  : "
                << time_.relativePath(outputName) << endl;
        }
    }

    return true;
}


// ************************************************************************* //
