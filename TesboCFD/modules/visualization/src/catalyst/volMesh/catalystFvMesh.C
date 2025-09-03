/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
    Copyright (C) 2018 CINECA
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

#include "catalystFvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{
    defineTypeNameAndDebug(fvMeshInput, 0);
    addNamedToRunTimeSelectionTable
    (
        catalystInput,
        fvMeshInput,
        dictionary,
        default
    );
}
} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::catalyst::fvMeshInput::update()
{
    // Backend requires a corresponding mesh
    backends_.filterKeys
    (
        [this](const word& k){ return meshes_.contains(k); }
    );

    forAllConstIters(meshes_, iter)
    {
        const word& regionName = iter.key();

        if (!backends_.contains(regionName))
        {
            auto backend =
                autoPtr<Foam::vtk::fvMeshAdaptor>::New
                (
                    *(iter.val()),
                    channelOpt_,
                    selectPatches_
                );

            // Special polyhedral treatment?
            backend->setDecompose(decomposeOpt_);

            backends_.set(regionName, backend);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyst::fvMeshInput::fvMeshInput
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    catalystInput(name),
    time_(runTime),
    channelOpt_(channelType::DEFAULT),
    decomposeOpt_(false),
    selectRegions_(),
    selectPatches_(),
    selectFields_(),
    meshes_(),
    backends_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalyst::fvMeshInput::read(const dictionary& dict)
{
    catalystInput::read(dict);

    meshes_.clear();
    backends_.clear();
    selectRegions_.clear();
    selectPatches_.clear();
    selectFields_.clear();
    decomposeOpt_ = dict.getOrDefault("decompose", false);

    unsigned selected(channelType::NONE);

    if (dict.getOrDefault("internal", true))
    {
        selected |= channelType::INTERNAL;
    }
    if (dict.getOrDefault("boundary", true))
    {
        selected |= channelType::BOUNDARY;
    }

    channelOpt_ = channelType(selected);

    if (channelType::NONE == channelOpt_)
    {
        // Both are disabled - warn, and skip the rest of the input.
        WarningInFunction
            << name() << ":" << nl
            << "    both internal and boundary are disabled" << nl
            << "    check your input" << nl;
        return true;
    }

    // All possible meshes
    meshes_ = time_.lookupClass<fvMesh>();

    dict.readIfPresent("patches", selectPatches_);
    dict.readIfPresent("regions", selectRegions_);

    if (selectRegions_.empty())
    {
        selectRegions_.resize(1);
        selectRegions_.front() =
            dict.getOrDefault<word>("region", polyMesh::defaultRegion);
    }

    // Restrict to specified meshes
    meshes_.filterKeys(selectRegions_);

    dict.readEntry("fields", selectFields_);

    return true;
}


void Foam::catalyst::fvMeshInput::update(polyMesh::readUpdateState state)
{
    // Trigger change of state

    // Be really paranoid and verify if the mesh actually exists
    const wordList regionNames(backends_.toc());

    for (const word& regionName : regionNames)
    {
        if (meshes_.contains(regionName) && time_.contains(regionName))
        {
            backends_[regionName]->updateState(state);
        }
        else
        {
            backends_.erase(regionName);
            meshes_.erase(regionName);
        }
    }
}


Foam::label Foam::catalyst::fvMeshInput::addChannels(dataQuery& dataq)
{
    if (channelType::NONE == channelOpt_)
    {
        // Everything disabled - do nothing
        return 0;
    }

    update();   // Enforce sanity for backends and adaptor

    if (backends_.empty())
    {
        return 0;
    }

    // Gather all fields that we know how to convert
    wordHashSet allFields;
    forAllConstIters(backends_, iter)
    {
        allFields += iter.val()->knownFields(selectFields_);
    }

    dataq.set(name(), allFields);

    return 1;
}


bool Foam::catalyst::fvMeshInput::convert
(
    dataQuery& dataq,
    outputChannels& outputs
)
{
    const word channelName(name());
    const wordList regionNames(backends_.sortedToc());

    if (regionNames.empty() || !dataq.contains(channelName))
    {
        // Not available, or not requested
        return false;
    }
    if (channelType::NONE == channelOpt_)
    {
        // Safety: everything disabled (should have been caught before)
        return false;
    }


    // A separate block for each region
    unsigned int blockNo = 0;

    for (const word& regionName : regionNames)
    {
        auto dataset = backends_[regionName]->output(selectFields_);

        // Existing or new
        vtkSmartPointer<vtkMultiBlockDataSet> block =
            outputs.lookup
            (
                channelName,
                vtkSmartPointer<vtkMultiBlockDataSet>::New()
            );

        block->SetBlock(blockNo, dataset);

        block->GetMetaData(blockNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            regionName                      // block name = region name
        );

        outputs.set(channelName, block);    // overwrites existing

        ++blockNo;
    }

    return true;
}


Foam::Ostream& Foam::catalyst::fvMeshInput::print(Ostream& os) const
{
    os  << name() << nl
        <<"    regions " << flatOutput(selectRegions_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl
        <<"    output ";

    const unsigned selected(channelOpt_);

    if (!selected)
    {
        os << " none";
    }
    else
    {
        if (0 != (selected & (channelType::INTERNAL)))
        {
            os  <<" internal";
        }
        if (0 != (selected & (channelType::BOUNDARY)))
        {
            os  <<" boundary";
        }
    }
    os << nl;

    return os;
}


// ************************************************************************* //
