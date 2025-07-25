/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "caseInfo.H"
#include "OFstream.H"
#include "fvMesh.H"
#include "cloud.H"
#include "globalMeshData.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvPatch.H"
#include "JSONformatter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(caseInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        caseInfo,
        dictionary
    );
}

const Enum<functionObjects::caseInfo::writeFormat>
Foam::functionObjects::caseInfo::writeFormatNames_
{
    { writeFormat::dict, "dictionary" },
    { writeFormat::json, "json" },
};

const Enum<functionObjects::caseInfo::lookupMode>
Foam::functionObjects::caseInfo::lookupModeNames_
{
    { lookupMode::none, "none" },
    { lookupMode::warn, "warn" },
    { lookupMode::error, "error" },
};
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::caseInfo::report(const string& str) const
{
    switch (lookupMode_)
    {
        case lookupMode::warn:
        {
            Warning << str.c_str() << endl;
            break;
        }
        case lookupMode::error:
        {
            FatalError << str.c_str() << exit(FatalError);
            break;
        }
        case lookupMode::none:
        {
            break;
        }
    }
}


void Foam::functionObjects::caseInfo::processDict
(
    dictionary& dict,
    const dictionary& targetDict,
    const entry* includePtr,
    const entry* excludePtr
) const
{
    auto sanitise = [](const wordRe& w){
        string str = w;
        str.replaceAll("/", "_").replaceAll(".", "_");

        // Strip any leading "_"
        while (str.starts_with('_'))
        {
            str = str.substr(1);
        }
        return str;
    };

    if (includePtr)
    {
        const wordRes includeEntryNames(includePtr->stream());
        for (const auto& nameRegex : includeEntryNames)
        {
            const auto* e = targetDict.findScoped(nameRegex, keyType::REGEX);
            if (e)
            {
                if (nameRegex.contains('/') || nameRegex.contains('.'))
                {
                    auto copyPtr = e->clone();
                    copyPtr->keyword() = sanitise(nameRegex);
                    dict.add(copyPtr.ptr());
                }
                else
                {
                    dict.add(*e);
                }
            }
            else
            {
                report
                (
                    "Unable to find entry "
                  + nameRegex
                  + " in dictionary "
                  + targetDict.name()
                );
            }
        }
    }
    else
    {
        if (excludePtr)
        {
            dictionary allData(targetDict);

            const wordRes excludeEntryNames(excludePtr->stream());

            for (const auto& nameRegex : excludeEntryNames)
            {
                const auto* e = allData.findScoped(nameRegex, keyType::REGEX);
                if (e)
                {
                    allData.remove(e->keyword());
                }
            }

            dict += allData;
        }
        else
        {
            // Add complete dictionary
            dict += targetDict;
        }
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::caseInfo::writeMeta(dictionary& out) const
{
    out.add("case", time_.globalCaseName());
    out.add("path", time_.globalPath());
    out.add("regions", time_.sortedNames<fvMesh>());
    out.add("nTimes", time_.times().size());
    out.add("nProc", Pstream::nProcs());
}


void Foam::functionObjects::caseInfo::writeRegisteredDicts
(
    const objectRegistry& obr,
    dictionary& out,
    dictionary& dictionaries
) const
{
    for (const auto& e : dictionaries)
    {
        const auto& keyword = e.keyword();

        if (!e.isDict())
        {
            FatalIOErrorInFunction(dictionaries)
                << "Entries must be specified in dictionary format. Please "
                << "correct entry " << keyword
                << exit(FatalIOError);
        }

        const dictionary& inputDict = e.dict();

        auto* includePtr = inputDict.findEntry("include");
        auto* excludePtr = inputDict.findEntry("exclude");

        const auto* ePtr = inputDict.findEntry("name");

        if (ePtr)
        {
            const word name(ePtr->stream());
            auto* dictPtr = obr.cfindObject<IOdictionary>(name);

            if (dictPtr)
            {
                processDict
                (
                    out.subDictOrAdd(keyword),
                   *dictPtr,
                    includePtr,
                    excludePtr
                );
                dictionaries.remove(keyword);
            }
        }
    }
}


void Foam::functionObjects::caseInfo::writeFileDicts
(
    dictionary& out,
    dictionary& dictionaries
) const
{
    for (auto& e : dictionaries)
    {
        const auto& keyword = e.keyword();

        if (!e.isDict())
        {
            FatalIOErrorInFunction(dictionaries)
                << "Entries must be specified in dictionary format. Please "
                << "correct entry " << keyword
                << exit(FatalIOError);
        }

        const dictionary& inputDict = e.dict();

        auto* includePtr = inputDict.findEntry("include");
        auto* excludePtr = inputDict.findEntry("exclude");

        const auto* ePtr = inputDict.findEntry("path");

        if (ePtr)
        {
            fileName path(ePtr->stream());
            path.expand();

            IOobject io(path, time_, IOobject::MUST_READ);

            if (!io.typeHeaderOk<dictionary>(false))
            {
                continue;
            }

            const word oldTypeName = IOdictionary::typeName;
            const_cast<word&>(IOdictionary::typeName) = word::null;

            processDict
            (
                out.subDictOrAdd(keyword),
                IOdictionary(io),
                includePtr,
                excludePtr
            );

            const_cast<word&>(IOdictionary::typeName) = oldTypeName;

            dictionaries.remove(keyword);
        }
    }
}


void Foam::functionObjects::caseInfo::writeFunctionObjects
(
    dictionary& out
) const
{
    for (const auto& fo : functionObjectNames_)
    {
        dictionary dict;
        if (getObjectResultDict(fo, dict))
        {
            out.add(fo, dict);
        }
        else
        {
            report("No result entries found for function object " + fo);
        }
    }
}


void Foam::functionObjects::caseInfo::writeMeshStats
(
    const polyMesh& mesh,
    dictionary& dict
) const
{
    dict.add("nGeometricD", mesh.nGeometricD());
    dict.add("nSolutionD", mesh.nSolutionD());

    const auto& globalData = mesh.globalData();

    dict.add("nPoints", globalData.nTotalPoints());
    dict.add("nFaces", globalData.nTotalFaces());
    dict.add("nCells", globalData.nTotalCells());

    dict.add("nPatches", mesh.boundaryMesh().nNonProcessor());

    dict.add("pointZones", mesh.pointZones().names());
    dict.add("faceZones", mesh.faceZones().names());
    dict.add("cellZones", mesh.cellZones().names());

    dict.add("boundsMin", mesh.bounds().min());
    dict.add("boundsMax", mesh.bounds().max());

    dict.add("clouds", mesh.sortedNames<cloud>());
}


namespace Foam
{
    template<class GeoFieldType>
    void addPatchTypeDetails(const fvMesh& mesh, dictionary& dict)
    {
        for (const GeoFieldType& obj : mesh.csorted<GeoFieldType>())
        {
            if (obj.readOpt() == IOobject::MUST_READ)
            {
                const auto& bf = obj.boundaryField();
                dictionary& objDict = dict.subDictOrAdd(obj.name());

                for (const auto& pf : bf)
                {
                    if (!isA<processorFvPatch>(pf.patch()))
                    {
                        objDict.add(pf.patch().name(), pf.type());
                    }
                }
            }
        }
    }

    template<template<typename> class FieldType>
    void addPatchDetails(const fvMesh& mesh, dictionary& dict)
    {
        addPatchTypeDetails<FieldType<scalar>>(mesh, dict);
        addPatchTypeDetails<FieldType<vector>>(mesh, dict);
        addPatchTypeDetails<FieldType<sphericalTensor>>(mesh, dict);
        addPatchTypeDetails<FieldType<symmTensor>>(mesh, dict);
        addPatchTypeDetails<FieldType<tensor>>(mesh, dict);
    }
}


void Foam::functionObjects::caseInfo::writePatches
(
    const fvMesh& mesh,
    dictionary& dict
) const
{
    // Geometry
    dictionary& bnd = dict.subDictOrAdd("types");
    const auto& pbm = mesh.boundaryMesh();
    for (const auto& pp : pbm)
    {
        if (!isA<processorPolyPatch>(pp))
        {
            bnd.add(pp.name(), pp.type());
        }
    }

    // Fields
    dictionary& fld = dict.subDictOrAdd("fields");
    addPatchDetails<VolumeField>(mesh, fld);
    addPatchDetails<SurfaceField>(mesh, fld);
    //addPatchDetails<AreaField>(mesh, fld);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::caseInfo::caseInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            name,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    stateFunctionObject(name, runTime),
    writeFile(runTime, name, typeName, dict),
    writeFormat_(writeFormat::dict),
    lookupMode_(lookupMode::warn)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::caseInfo::read(const dictionary& dict)
{
    if (stateFunctionObject::read(dict) && writeFile::read(dict))
    {
        writeFormatNames_.readIfPresent("writeFormat", dict, writeFormat_);
        writeFile::setExt("." + writeFormatNames_[writeFormat_]);

        lookupModeNames_.readIfPresent("lookupMode", dict, lookupMode_);

        dictionaries_ = dict.subOrEmptyDict("dictionaries");

        dict.readIfPresent("functionObjects", functionObjectNames_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::caseInfo::execute()
{
    return true;
}


bool Foam::functionObjects::caseInfo::write()
{
    // Output dictionary
    dictionary& data = *this;
    data.clear();

    // Case meta data
    writeMeta(data.subDictOrAdd("meta"));

    // Note: copying dictionaries
    // - these are consumed/removed when found to enable checks that all
    //   dictionaries are processed
    dictionary dicts(dictionaries_);

    dictionary& dataDicts = data.subDictOrAdd("dictionaries");

    // Time-registered dictionaries
    writeRegisteredDicts(time_, dataDicts, dicts);

    // File-based dictionaries
    writeFileDicts(dataDicts, dicts);

    // Per-region information
    dictionary& regionDict = data.subDictOrAdd("regions");

    for (const fvMesh& mesh : time_.csorted<fvMesh>())
    {
        dictionary meshDicts(dicts);

        const word& name = mesh.name();

        dictionary& out = regionDict.subDictOrAdd(name);

        writeMeshStats(mesh, out.subDictOrAdd("mesh"));

        writePatches(mesh, out.subDictOrAdd("boundary"));

        // Mesh-registered dictionaries
        writeRegisteredDicts(mesh, out.subDictOrAdd("dictionaries"), meshDicts);

        for (const word& keyword : meshDicts.sortedToc())
        {
            report
            (
                "Mesh '"
               + keyword
               + "' : Unable to process dictionary entry '"
               + keyword
               + "'"
            );
        }
    }


    writeFunctionObjects(data.subDictOrAdd("functions"));


    if (Pstream::master())
    {
        auto filePtr = newFileAtTime(functionObject::name(), time_.value());
        auto& os = filePtr();

        // Reset stream width - was set in writeFile
        os.width(0);

        switch (writeFormat_)
        {
            case writeFormat::dict:
            {
                IOdictionary::writeData(os);
                break;
            }
            case writeFormat::json:
            {
                JSONformatter json(os);
                json.writeDict(data);
                break;
            }
        }

        Info<< "Written " << writeFormatNames_[writeFormat_]
            << " file: " << os.name() << endl;
    }

    return true;
}


// ************************************************************************* //
