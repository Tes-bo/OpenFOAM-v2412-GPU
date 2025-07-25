/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "ensightWrite.H"
#include "ensightOutput.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ensightWrite, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ensightWrite,
        dictionary
    );
}
}

// Implementation
#include "ensightWriteImpl.C"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::ensightWrite::writeAllVolFields
(
    const fvMeshSubset& proxy,
    const wordHashSet& candidateNames
)
{
    label count = 0;

    ensightOutput::floatBufferType scratch;

    {
        #undef  doLocalCode
        #define doLocalCode(PrimitiveType)              \
            count += writeVolFieldsImpl<PrimitiveType>  \
            (                                           \
                scratch,                                \
                proxy,                                  \
                candidateNames                          \
            );

        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }

    return count;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ensightWrite::ensightWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeOpts_(),
    caseOpts_("format", dict, IOstreamOption::BINARY),
    outputDir_(),
    consecutive_(false),
    meshState_(polyMesh::TOPO_CHANGE),
    selectFields_(),
    blockFields_(),
    selection_(),
    meshSubset_(mesh_),
    ensCase_(nullptr),
    ensMesh_(nullptr)
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

bool Foam::functionObjects::ensightWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    readSelection(dict);


    // Writer options

    consecutive_ = dict.getOrDefault("consecutive", false);

    writeOpts_.useBoundaryMesh(dict.getOrDefault("boundary", true));
    writeOpts_.useInternalMesh(dict.getOrDefault("internal", true));


    // Warn if noPatches keyword (1806) exists and contradicts our settings
    // Cannot readily use Compat since the boolean has the opposite value.
    if
    (
        dict.getOrDefault("noPatches", false)
     && writeOpts_.useBoundaryMesh()
    )
    {
        WarningInFunction
            << "Use 'boundary' instead of 'noPatches' to enable/disable "
            << "conversion of the boundaries" << endl;
    }

    wordRes list;
    if (dict.readIfPresent("patches", list))
    {
        list.uniq();
        writeOpts_.patchSelection(list);
    }
    if (dict.readIfPresent("excludePatches", list))
    {
        list.uniq();
        writeOpts_.patchExclude(list);
    }

    if (dict.readIfPresent("faceZones", list))
    {
        list.uniq();
        writeOpts_.faceZoneSelection(list);
    }


    // Case options

    caseOpts_.nodeValues(dict.getOrDefault("nodeValues", false));
    caseOpts_.width(dict.getOrDefault<label>("width", 8));
    caseOpts_.overwrite(dict.getOrDefault("overwrite", false));

    caseOpts_.timeFormat("timeFormat", dict);
    caseOpts_.timePrecision("timePrecision", dict);


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


bool Foam::functionObjects::ensightWrite::execute()
{
    return true;
}


bool Foam::functionObjects::ensightWrite::write()
{
    if (!ensCase_)
    {
        ensCase_.reset
        (
            new ensightCase(outputDir_, time_.globalCaseName(), caseOpts_)
        );
    }

    if (consecutive_)
    {
        ensCase().nextTime(time_.value());
    }
    else
    {
        ensCase().setTime(time_.value(), time_.timeIndex());
    }


    if (update())
    {
        // Treat all geometry as moving, since we do not know a priori
        // if the simulation has mesh motion later on.
        autoPtr<ensightGeoFile> os = ensCase_().newGeometry(true);
        ensMesh_().write(os.ref());
    }


    // Output fields MUST be specified to avoid accidentally
    // writing everything. Can still use ".*" for everything

    wordHashSet candidateNames;

    if (!selectFields_.empty())
    {
        if (!blockFields_.empty())
        {
            // With 'allow' and 'deny' filters
            wordRes::filter filter(selectFields_, blockFields_);

            candidateNames = mesh_.names<void>(filter);
        }
        else
        {
            // With 'allow' filter only
            candidateNames = mesh_.names<void>(selectFields_);
        }
    }

    // Prune restart fields
    candidateNames.filterKeys
    (
        [](const word& k){ return k.ends_with("_0"); },
        true // prune
    );

    Log << type() << " " << name() << " write:\n";
    writeAllVolFields(meshSubset_, candidateNames);

    Log << nl;

    ensCase().write();  // Flush case information

    return true;
}


bool Foam::functionObjects::ensightWrite::end()
{
    return true;
}


// ************************************************************************* //
