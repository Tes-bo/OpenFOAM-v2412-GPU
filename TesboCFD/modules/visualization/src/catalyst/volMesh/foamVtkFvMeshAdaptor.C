/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "foamVtkFvMeshAdaptor.H"
#include "Pstream.H"

// OpenFOAM includes
#include "fvMesh.H"

// VTK includes
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{
    defineTypeNameAndDebug(fvMeshAdaptor, 0);
}
} // End namespace Foam


const Foam::Enum
<
    Foam::vtk::fvMeshAdaptor::channelType
>
Foam::vtk::fvMeshAdaptor::channelNames
{
    { channelType::NONE,     "none" },
    { channelType::INTERNAL, "internal" },
    { channelType::BOUNDARY, "boundary" },
    { channelType::ALL,      "all" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::fvMeshAdaptor::definePatchIds()
{
    // Generate or update the list of patchIds

    patchIds_.clear();

    if (!usingBoundary())
    {
        return;
    }

    // General patch information
    // Restrict to non-processor patches.
    // This value is invariant across all processors.

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nNonProcessor = pbm.nNonProcessor();

    if (patchPatterns_.empty())
    {
        patchIds_ = identity(nNonProcessor);
    }
    else
    {
        // Don't warn if not found, use patch groups
        labelHashSet ids(pbm.patchSet(patchPatterns_, false, true));

        // Restricted to non-processor patches
        ids.filterKeys
        (
            [nNonProcessor](label i){ return i < nNonProcessor; }
        );

        // MUST be sorted. Other internal logic relies upon this!
        patchIds_ = ids.sortedToc();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::fvMeshAdaptor::fvMeshAdaptor
(
    const fvMesh& mesh,
    const channelType channelsOpt,
    const wordRes& patchSelection
)
:
    mesh_(mesh),
    patchPatterns_(patchSelection),
    patchIds_(),
    channels_(channelsOpt),
    interpFields_(true),
    extrapPatches_(false),
    decomposePoly_(false),
    meshState_(polyMesh::TOPO_CHANGE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::fvMeshAdaptor::setChannels(const wordList& chanNames)
{
    unsigned chanIds = 0;
    for (const word& chan : chanNames)
    {
        if (channelNames.contains(chan))
        {
            chanIds |= channelNames[chan];
        }
    }

    setChannels(chanIds);
}


void Foam::vtk::fvMeshAdaptor::setChannels(enum channelType chanIds)
{
    channels_ = chanIds;

    if (!usingInternal())
    {
        cachedVtu_.clear();
    }

    if (!usingBoundary())
    {
        cachedVtp_.clear();
        patchIds_.clear();
    }
}


void Foam::vtk::fvMeshAdaptor::setChannels(unsigned chanIds)
{
    channels_ = (chanIds & 0x3);

    if (!usingInternal())
    {
        cachedVtu_.clear();
    }

    if (!usingBoundary())
    {
        cachedVtp_.clear();
        patchIds_.clear();
    }
}


void Foam::vtk::fvMeshAdaptor::setDecompose(const bool val)
{
    if (usingInternal() && val != decomposePoly_)
    {
        cachedVtu_.clear();
        decomposePoly_ = val;
    }
}


Foam::label Foam::vtk::fvMeshAdaptor::channels() const
{
    return label(channels_);
}


bool Foam::vtk::fvMeshAdaptor::usingInternal() const
{
    return (channels_ & INTERNAL);
}


bool Foam::vtk::fvMeshAdaptor::usingBoundary() const
{
    return (channels_ & BOUNDARY);
}


const Foam::labelList& Foam::vtk::fvMeshAdaptor::patchIds() const
{
    return patchIds_;
}


void Foam::vtk::fvMeshAdaptor::updateContent(const wordRes& selectFields)
{
    const bool oldDecomp = decomposePoly_;

    // General patch information
    // Restrict to non-processor patches.
    // This value is invariant across all processors.

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nNonProcessor = pbm.nNonProcessor();


    // Update cached, saved, unneed values.

    HashSet<string> nowActive;

    // INTERNAL
    if (usingInternal())
    {
        nowActive.insert(internalName());
    }


    // BOUNDARY
    if (usingBoundary())
    {
        for (label patchId=0; patchId < nNonProcessor; ++patchId)
        {
            const polyPatch& pp = pbm[patchId];
            nowActive.insert(pp.name());
        }
    }

    // Dispose of unneeded components
    cachedVtp_.retain(nowActive);
    cachedVtu_.retain(nowActive);

    if
    (
        meshState_ == polyMesh::TOPO_CHANGE
     || meshState_ == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Eliminate cached values that would be unreliable
        forAllIters(cachedVtp_, iter)
        {
            iter.val().clearGeom();
            iter.val().clear();
        }
        forAllIters(cachedVtu_, iter)
        {
            iter.val().clearGeom();
            iter.val().clear();
        }

        definePatchIds();
    }
    else if (oldDecomp != decomposePoly_)
    {
        // poly-decompose changed - dispose of cached values
        forAllIters(cachedVtu_, iter)
        {
            iter.val().clearGeom();
            iter.val().clear();
        }
    }

    convertGeometryInternal();
    convertGeometryBoundary();

    applyGhosting();
    convertVolFields(selectFields);
    meshState_ = polyMesh::UNCHANGED;
}


vtkSmartPointer<vtkMultiBlockDataSet>
Foam::vtk::fvMeshAdaptor::output(const wordRes& select)
{
    updateContent(select);

    // All individual datasets are vtkMultiPieceDataSet for improved
    // handling downstream.

    label rank = 0;
    label nproc = 1;

    if (UPstream::parRun())
    {
        rank  = UPstream::myProcNo();
        nproc = UPstream::nProcs();
    }


    auto outputs = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    unsigned int blockNo = 0;

    // INTERNAL
    if (usingInternal())
    {
        do
        {
            const auto& longName = internalName();
            auto iter = cachedVtu_.find(longName);
            if (!iter.good() || !iter.val().dataset)
            {
                Pout<<"Cache miss for VTU " << longName << endl;
                break; // Should never happen
            }

            foamVtuData& vtuData = iter.val();

            auto pieces = vtkSmartPointer<vtkMultiPieceDataSet>::New();

            pieces->SetNumberOfPieces(nproc);
            pieces->SetPiece(rank, vtuData.dataset);

            outputs->SetBlock(blockNo, pieces);

            outputs->GetMetaData(blockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                internalName()
            );

            ++blockNo;
        } while (false);  // do once
    }

    // BOUNDARY
    if (!patchIds_.empty())
    {
        unsigned int subBlockNo = 0;

        auto output = vtkSmartPointer<vtkMultiBlockDataSet>::New();

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        for (const label patchId : patchIds_)
        {
            const polyPatch& pp = pbm[patchId];
            const word& longName = pp.name();

            auto iter = cachedVtp_.find(longName);
            if (!iter.good() || !iter.val().dataset)
            {
                Pout<<"Cache miss for VTP patch " << longName << endl;
                break; // Should never happen
            }

            foamVtpData& vtpData = iter.val();

            auto pieces = vtkSmartPointer<vtkMultiPieceDataSet>::New();

            pieces->SetNumberOfPieces(nproc);
            pieces->SetPiece(rank, vtpData.dataset);

            output->SetBlock(subBlockNo, pieces);

            output->GetMetaData(subBlockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                longName
            );

            ++subBlockNo;
        }

        outputs->SetBlock(blockNo, output);

        outputs->GetMetaData(blockNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            boundaryName()
        );
    }

    // Would actually like to have this:
    //     outputs->SetName(mesh_.name().c_str());
    // but do that in the caller side

    return outputs;
}


void Foam::vtk::fvMeshAdaptor::updateState(polyMesh::readUpdateState state)
{
    // Only move to worse states
    switch (state)
    {
        case polyMesh::UNCHANGED:
            break;

        case polyMesh::POINTS_MOVED:
            if (meshState_ == polyMesh::UNCHANGED)
            {
                meshState_ = polyMesh::POINTS_MOVED;
            }
            break;

        case polyMesh::TOPO_CHANGE:
        case polyMesh::TOPO_PATCH_CHANGE:
            meshState_ = polyMesh::TOPO_CHANGE;
            break;
    }
}


// ************************************************************************* //
