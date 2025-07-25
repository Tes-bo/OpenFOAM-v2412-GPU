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

#include "vtkCloud.H"
#include "Cloud.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "foamVtkOutputOptions.H"
#include "addToRunTimeSelectionTable.H"
#include "pointList.H"
#include "stringOps.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vtkCloud, 0);
    addToRunTimeSelectionTable(functionObject, vtkCloud, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::vtkCloud::writeVerts
(
    autoPtr<vtk::formatter>& format,
    const label nTotParcels
) const
{
    // No collectives - can skip on slave processors
    if (!format) return;

    // Same payload for connectivity and offsets
    const uint64_t payLoad = vtk::sizeofData<label>(nTotParcels);

    format().tag(vtk::fileTag::VERTS);

    //
    // 'connectivity'
    // = linear mapping onto points
    //
    {
        format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
        format().writeSize(payLoad);

        vtk::writeIdentity(format(), nTotParcels);

        format().flush();
        format().endDataArray();
    }

    //
    // 'offsets' (connectivity offsets)
    // = linear mapping onto points (with 1 offset)
    //
    {
        format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
        format().writeSize(payLoad);

        vtk::writeIdentity(format(), nTotParcels, 1);

        format().flush();
        format().endDataArray();
    }

    format().endTag(vtk::fileTag::VERTS);
}


bool Foam::functionObjects::vtkCloud::writeCloud
(
    const fileName& file,
    const word& cloudName
)
{
    applyFilter_ = false;

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
            "vtk::vtkCloud::" + cloudName,
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

    if (pruneEmpty_ && !nTotParcels)
    {
        return false;
    }

    std::ofstream os;
    autoPtr<vtk::formatter> format;

    if (!file.has_ext("vtp"))
    {
        FatalErrorInFunction
            << type() << " File missing .vtp extension!" << nl << endl
            << exit(FatalError);
    }

    if (UPstream::master())
    {
        Foam::mkDir(file.path());
        os.open(file);

        format = writeOpts_.newFormatter(os);

        // beginFile()

        // XML (inline)
        format()
            .xmlHeader()
            .xmlComment
            (
                "case='" + time_.globalCaseName()
              + "' cloud='" + cloudName
              + "' time='" + time_.timeName()
              + "' index='" + Foam::name(time_.timeIndex())
              + "'"
            )
            .beginVTKFile<vtk::fileTag::POLY_DATA>();


        // FieldData with TimeValue
        format()
            .beginFieldData()
            .writeTimeValue(time_.value())
            .endFieldData();


        // writeGeometry()

        // beginPiece()
        if (useVerts_)
        {
            format()
                .tag
                (
                    vtk::fileTag::PIECE,
                    vtk::fileAttr::NUMBER_OF_POINTS, nTotParcels,
                    vtk::fileAttr::NUMBER_OF_VERTS, nTotParcels
                );
        }
        else
        {
            format()
                .tag
                (
                    vtk::fileTag::PIECE,
                    vtk::fileAttr::NUMBER_OF_POINTS, nTotParcels
                );
        }

        // writePoints()
        {
            const uint64_t payLoad = vtk::sizeofData<float,3>(nTotParcels);

            format().tag(vtk::fileTag::POINTS)
                .beginDataArray<float,3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
        }
    }


    if (applyFilter_)
    {
        vtk::writeListParallel(format.ref(), *pointsPtr, parcelAddr_);
    }
    else
    {
        vtk::writeListParallel(format.ref(), *pointsPtr);
    }


    if (UPstream::master())
    {
        format().flush();
        format().endDataArray();
        format().endTag(vtk::fileTag::POINTS);

        if (useVerts_)
        {
            writeVerts(format, nTotParcels);
        }
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

    // Restrict to specified fields
    if (selectFields_.size())
    {
        obrTmp.filterKeys(selectFields_);
    }


    // Write fields

    if (UPstream::master())
    {
        if (useVerts_)
        {
            format().beginCellData();
        }
        else
        {
            format().beginPointData();
        }
    }

    DynamicList<word> written(obrTmp.size() + currCloud.objectRegistry::size());

    written.push_back
    (
        writeFields<label>(format, obrTmp, nTotParcels)
    );
    written.push_back
    (
        writeFields<scalar>(format, obrTmp, nTotParcels)
    );
    written.push_back
    (
        writeFields<vector>(format, obrTmp, nTotParcels)
    );

    // Any cloudFunctions results
    written.push_back
    (
        writeFields<scalar>(format, currCloud, nTotParcels)
    );

    if (UPstream::master())
    {
        if (useVerts_)
        {
            format().endCellData();
        }
        else
        {
            format().endPointData();
        }

        format().endPiece();
        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }


    // Record information into the state (all processors)
    //
    // foName
    // {
    //     cloudName
    //     {
    //         file   "<case>/postProcessing/name/cloud1_0001.vtp";
    //         fields (U T rho);
    //     }
    // }

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

Foam::functionObjects::vtkCloud::vtkCloud
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeOpts_(vtk::formatType::INLINE_BASE64),
    printf_(),
    useVerts_(false),
    pruneEmpty_(false),
    applyFilter_(false),
    selectClouds_(),
    selectFields_(),
    directory_(),
    series_()
{
    // May still want this? (OCT-2018)
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

bool Foam::functionObjects::vtkCloud::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // We probably cannot trust old information after a reread
    series_.clear();

    //
    // Default format is xml base64. Legacy is not desired.
    //
    writeOpts_ = vtk::formatType::INLINE_BASE64;

    writeOpts_.ascii
    (
        IOstreamOption::ASCII
     == IOstreamOption::formatEnum("format", dict, IOstreamOption::BINARY)
    );

    writeOpts_.append(false);  // No append supported
    writeOpts_.legacy(false);  // No legacy supported

    writeOpts_.precision
    (
        dict.getOrDefault("precision", IOstream::defaultPrecision())
    );

    // Info<< type() << " " << name() << " output-format: "
    //     << writeOpts_.description() << nl;

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

    // useTimeName_ = dict.getOrDefault("useTimeName", false);

    useVerts_ = dict.getOrDefault("cellData", false);
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


bool Foam::functionObjects::vtkCloud::execute()
{
    return true;
}


bool Foam::functionObjects::vtkCloud::write()
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

    const scalar timeValue = time_.value();

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
            directory_/cloudName + timeDesc + ".vtp"
        );

        // writeCloud() includes mkDir (on master)

        if (writeCloud(outputName, cloudName))
        {
            Log << "    cloud  : "
                << time_.relativePath(outputName) << endl;

            if (UPstream::master())
            {
                // Add to file-series and emit as JSON
                fileName seriesName(vtk::seriesWriter::base(outputName));

                vtk::seriesWriter& series = series_(seriesName);

                // First time?
                // Load from file, verify against filesystem,
                // prune time >= currentTime
                if (series.empty())
                {
                    series.load(seriesName, true, timeValue);
                }

                series.append(timeValue, outputName);
                series.write(seriesName);
            }
        }
    }

    return true;
}


// ************************************************************************* //
