/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "catalystFaMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"
#include "fvMesh.H"

#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{
    defineTypeNameAndDebug(faMeshInput, 0);
    addNamedToRunTimeSelectionTable
    (
        catalystInput,
        faMeshInput,
        dictionary,
        area
    );
}
} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::catalyst::faMeshInput::update()
{
    // Backend requires a corresponding mesh
    backends_.filterKeys
    (
        [this](const word& k){ return meshes_.contains(k); }
    );

    forAllConstIters(meshes_, iter)
    {
        const word& name = iter.key();
        const auto* mesh = iter.val();

        if (!backends_.contains(name))
        {
            auto backend = autoPtr<Foam::vtk::faMeshAdaptor>::New(*mesh);

            // Apply any configuration options
            // ...

            backends_.set(name, backend);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyst::faMeshInput::faMeshInput
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    catalystInput(name),
    time_(runTime),
    regionName_(),
    selectAreas_(),
    selectFields_(),
    meshes_(),
    backends_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalyst::faMeshInput::read(const dictionary& dict)
{
    catalystInput::read(dict);

    selectAreas_.clear();
    selectFields_.clear();
    backends_.clear();

    regionName_ = dict.getOrDefault<word>("region", polyMesh::defaultRegion);

    const polyMesh& pMesh = time_.lookupObject<polyMesh>(regionName_);

    // Registry containing all finite-area meshes on the polyMesh
    const auto* faRegistry = faMesh::registry(pMesh);

    if (faRegistry)
    {
        // All possible meshes for the given region
        meshes_ = faRegistry->lookupClass<faMesh>();
    }
    else
    {
        meshes_.clear();
    }

    dict.readIfPresent("areas", selectAreas_);

    if (selectAreas_.empty())
    {
        word areaName;
        if (!dict.readIfPresent("area", areaName))
        {
            wordList available = meshes_.sortedToc();

            if (!available.empty())
            {
                areaName = available.front();
            }
        }

        if (!areaName.empty())
        {
            selectAreas_.resize(1);
            selectAreas_.front() = areaName;
        }
    }

    // Restrict to specified meshes
    meshes_.filterKeys(selectAreas_);

    dict.readEntry("fields", selectFields_);

    return true;
}


void Foam::catalyst::faMeshInput::update(polyMesh::readUpdateState state)
{
    // Trigger change of state

    const polyMesh& pMesh = time_.lookupObject<polyMesh>(regionName_);

    // Registry containing all finite-area meshes on the polyMesh
    const auto* faRegistry = faMesh::registry(pMesh);

    // Be really paranoid and verify if the mesh actually exists
    const wordList areaNames(backends_.sortedToc());

    for (const word& areaName : areaNames)
    {
        if
        (
            meshes_.contains(areaName)
         && faRegistry && faRegistry->cfindObject<faMesh>(areaName)
        )
        {
            backends_[areaName]->updateState(state);
        }
        else
        {
            backends_.erase(areaName);
            meshes_.erase(areaName);
        }
    }
}


Foam::label Foam::catalyst::faMeshInput::addChannels(dataQuery& dataq)
{
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


bool Foam::catalyst::faMeshInput::convert
(
    dataQuery& dataq,
    outputChannels& outputs
)
{
    const word channelName(name());
    const wordList areaNames(backends_.sortedToc());

    if (areaNames.empty() || !dataq.contains(channelName))
    {
        // Not available, or not requested
        return false;
    }


    // A separate block for each area mesh
    unsigned int blockNo = 0;

    for (const word& areaName : areaNames)
    {
        auto dataset = backends_[areaName]->output(selectFields_);

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
            areaName                        // block name = area mesh name
        );

        outputs.set(channelName, block);    // overwrites existing

        ++blockNo;
    }

    return true;
}


Foam::Ostream& Foam::catalyst::faMeshInput::print(Ostream& os) const
{
    os  << name() << nl
        <<"    areas   " << flatOutput(selectAreas_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl;

    return os;
}


// ************************************************************************* //
