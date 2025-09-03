/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "solverTemplate.H"
#include "Time.H"
#include "IOPtrList.H"
#include "polyMesh.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::solverTemplate::solverType
>
Foam::solverTemplate::solverTypeNames_
({
    { solverType::stCompressible, "compressible" },
    { solverType::stIncompressible, "incompressible" },
    { solverType::stBuoyant, "buoyant" },
    { solverType::stUnknown, "unknown" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::solverTemplate::readFromDict
(
    IOobject& dictHeader,
    const word& entryName
) const
{
    if (!dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << "Unable to open file "
            << dictHeader.objectPath()
            << exit(FatalError);
    }

    IOdictionary dict(dictHeader);
    return dict.get<word>(entryName);
}


Foam::dictionary Foam::solverTemplate::readFluidFieldTemplates
(
    const word& regionName,
    const fileName& baseDir,
    const dictionary& solverDict,
    const Time& runTime
) const
{
    Info<< "    Reading fluid field templates";

    if (!regionName.empty())
    {
        Info<< " for region " << regionName;
    }
    Info<< endl;

    dictionary fieldTemplates = solverDict.subDict("fluidFields");

    const fileName turbModelDir(baseDir/"models"/"turbulence");

    const dictionary fieldModels(solverDict.subDict("fluidModels"));

    word turbModel("laminar"); // default to laminar
    word turbType("none");

    if (fieldModels.readIfPresent("turbulenceModel", turbType))
    {
        if (turbType == "turbulenceModel")
        {
            IOdictionary turbPropDict
            (
                IOobject
                (
                    // turbulenceModel::propertiesName
                    "turbulenceProperties",
                    runTime.constant(),
                    regionName,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            const word modelType(turbPropDict.get<word>("simulationType"));

            if (modelType == "laminar")
            {
                // Leave as laminar
            }
            else if (modelType == "RAS")
            {
                // "RASModel" for v2006 and earlier
                turbPropDict.subDict(modelType)
                    .readCompat("model", {{"RASModel", -2006}}, turbModel);
            }
            else if (modelType == "LES")
            {
                // "LESModel" for v2006 and earlier
                turbPropDict.subDict(modelType)
                    .readCompat("model", {{"LESModel", -2006}}, turbModel);
            }
            else
            {
                FatalErrorInFunction
                    << "Unhandled turbulence model option " << modelType
                    << ". Valid options are laminar, RAS, LES"
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unhandled turbulence model option " << turbType
                << ". Valid options are turbulenceModel"
                << exit(FatalError);
        }
    }

    Info<< "    Selecting " << turbType << ": " << turbModel << endl;

    IOdictionary turbModelDict
    (
        IOobject
        (
            fileName(turbModelDir/turbModel),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Merge common fluid fields
    fieldTemplates.merge(turbModelDict.subDict("fluidFields"));

    // Merge specific compressible or incompressible fluid fields
    switch (solverType_)
    {
        case stIncompressible:
        {
            fieldTemplates.merge(turbModelDict.subDict("incompressibleFields"));
            break;
        }
        case stCompressible:
        case stBuoyant:
        {
            fieldTemplates.merge(turbModelDict.subDict("compressibleFields"));
            break;
        }
        default:
        {
            // do nothing
        }
    }

    return fieldTemplates;
}


Foam::dictionary Foam::solverTemplate::readSolidFieldTemplates
(
    const word& regionName,
    const dictionary& solverDict
) const
{
    Info<< "    Reading solid field templates for region " << regionName
        << endl;

    // nothing more to do for solids
    return solverDict.subDict("solidFields");
}


void Foam::solverTemplate::setRegionProperties
(
    const dictionary& fieldDict,
    const word& regionType,
    const word& regionName,
    const label regionI
)
{
    regionTypes_[regionI] = regionType,
    regionNames_[regionI] = regionName,
    fieldNames_[regionI] = fieldDict.toc();
    fieldTypes_[regionI].setSize(fieldNames_[regionI].size());
    fieldDimensions_[regionI].setSize(fieldNames_[regionI].size());

    forAll(fieldNames_[regionI], i)
    {
        const word& fieldName = fieldNames_[regionI][i];
        const dictionary& dict = fieldDict.subDict(fieldName);

        dict.readEntry("type", fieldTypes_[regionI][i]);
        fieldDimensions_[regionI].set
        (
            i,
            new dimensionSet("dimensions", dict)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solverTemplate::solverTemplate
(
    const fileName& baseDir,
    const Time& runTime,
    const word& solverName
)
:
    solverType_(stUnknown),
    multiRegion_(false),
    regionTypes_(),
    fieldNames_(),
    fieldTypes_(),
    fieldDimensions_()
{
    IOdictionary solverDict
    (
        IOobject
        (
            fileName(baseDir/"solvers"/solverName),
            runTime,
            IOobject::MUST_READ
        )
    );

    Info<< "Selecting " << solverName << ": ";

    solverType_ = solverTypeNames_.get("solverType", solverDict);
    multiRegion_ = solverDict.get<bool>("multiRegion");

    Info<< solverTypeNames_[solverType_];
    if (multiRegion_)
    {
        Info<< ", multi-region";
    }
    else
    {
        Info<< ", single-region";
    }

    Info<< " case" << endl;

    if (multiRegion_)
    {
        // read regionProperties
        regionProperties rp(runTime);

        const wordList fluidNames(rp["fluid"]);
        const wordList solidNames(rp["solid"]);

        const label nRegion = fluidNames.size() + solidNames.size();

        regionTypes_.setSize(nRegion);
        regionNames_.setSize(nRegion);
        fieldNames_.setSize(nRegion);
        fieldTypes_.setSize(nRegion);
        fieldDimensions_.setSize(nRegion);

        // read templates for solver lists of available
        // - fields and dimensions required for the solver
        // - models
        label regionI = 0;
        forAll(fluidNames, i)
        {
            const dictionary fieldDict
            (
                readFluidFieldTemplates
                (
                    fluidNames[i],
                    baseDir,
                    solverDict,
                    runTime
                )
            );

            setRegionProperties(fieldDict, "fluid", fluidNames[i], regionI++);
        }

        forAll(solidNames, i)
        {
            const dictionary fieldDict
            (
                readSolidFieldTemplates
                (
                    solidNames[i],
                    solverDict
                )
            );
            setRegionProperties(fieldDict, "solid", solidNames[i], regionI++);
        }
    }
    else
    {
        regionTypes_.setSize(1);
        regionNames_.setSize(1);
        fieldNames_.setSize(1);
        fieldTypes_.setSize(1);
        fieldDimensions_.setSize(1);

        // read templates for solver lists of available
        // - fields and dimensions required for the solver
        // - models
        const dictionary fieldDict
        (
            readFluidFieldTemplates
            (
                word::null, // use region = null for single-region cases
                baseDir,
                solverDict,
                runTime
            )
        );

        setRegionProperties(fieldDict, "fluid", word::null, 0);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::solverTemplate::type() const
{
    return solverTypeNames_[solverType_];
}


bool Foam::solverTemplate::multiRegion() const
{
    return multiRegion_;
}


Foam::label Foam::solverTemplate::nRegion() const
{
    return regionTypes_.size();
}


const Foam::word& Foam::solverTemplate::regionType(const label regionI) const
{
    return regionTypes_[regionI];
}


const Foam::word& Foam::solverTemplate::regionName(const label regionI) const
{
    return regionNames_[regionI];
}


const Foam::wordList& Foam::solverTemplate::fieldNames
(
    const label regionI
) const
{
    return fieldNames_[regionI];
}


const Foam::wordList& Foam::solverTemplate::fieldTypes
(
    const label regionI
) const
{
    return fieldTypes_[regionI];
}


const Foam::PtrList<Foam::dimensionSet>& Foam::solverTemplate::fieldDimensions
(
    const label regionI
) const
{
    return fieldDimensions_[regionI];
}


// ************************************************************************* //
