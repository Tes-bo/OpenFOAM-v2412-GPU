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
#include "cellCellStencilObject.H"

// VTK includes
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::fvMeshAdaptor::convertGeometryInternal()
{
    // INTERNAL
    if (!usingInternal())
    {
        cachedVtu_.clear();
        return;
    }

    const auto& longName = internalName();

    foamVtuData& vtuData = cachedVtu_(longName);

    vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
    if (vtuData.nPoints())
    {
        if (meshState_ == polyMesh::UNCHANGED)
        {
            DebugInfo
                << "reuse " << longName << nl;

            vtuData.reuse();  // No movement - simply reuse
            return;
        }
        else if (meshState_ == polyMesh::POINTS_MOVED)
        {
            DebugInfo
                << "move points " << longName << nl;

            vtkgeom = vtuData.getCopy();
            vtkgeom->SetPoints(vtuData.points(mesh_));
        }
    }

    if (!vtkgeom)
    {
        DebugInfo
            << "Nothing usable from cache - create new geometry" << nl;

        // Nothing usable from cache - create new geometry
        vtkgeom = vtuData.internal(mesh_, decomposePoly_);
    }

    vtuData.set(vtkgeom);
}


void Foam::vtk::fvMeshAdaptor::convertGeometryBoundary()
{
    // BOUNDARY

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const label patchId : patchIds_)
    {
        const polyPatch& pp = patches[patchId];
        const word& longName = pp.name();

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo
                    << "reuse " << longName << nl;

                vtpData.reuse();  // No movement - simply reuse
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Point movement on single patch
                DebugInfo
                    << "move points " << longName << nl;

                vtkgeom = vtpData.getCopy();
                vtkgeom->SetPoints(vtk::Tools::Patch::points(pp));
            }
        }

        vtpData.clear(); // Remove any old mappings

        DebugInfo
            << "Creating VTK mesh for patch [" << patchId <<"] "
            << longName << endl;

        // Unused information
        vtpData.additionalIds().clear();

        // This is somewhat inconsistent, since we currently only have
        // normal (non-grouped) patches but this may change in the future.

        vtkgeom = vtk::Tools::Patch::mesh(pp);

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            // Catch any problems
            cachedVtp_.erase(longName);
        }
    }
}


void Foam::vtk::fvMeshAdaptor::applyGhostingInternal(const labelUList& types)
{
    // INTERNAL

    if (types.empty() || !usingInternal())
    {
        return;
    }

    const auto& longName = internalName();

    auto iter = cachedVtu_.find(longName);
    if (!iter.good() || !iter.val().dataset)
    {
        // Should not happen, but for safety require a vtk geometry
        Pout<<"Cache miss for VTU " << longName << endl;
        return;
    }
    foamVtuData& vtuData = iter.val();
    auto dataset = vtuData.dataset;

    const labelUList& cellMap = vtuData.cellMap();

    auto vtkgcell = dataset->GetCellGhostArray();
    if (!vtkgcell)
    {
        vtkgcell = dataset->AllocateCellGhostArray();
    }

    // auto vtkgpoint = dataset->GetPointGhostArray();
    // if (!vtkgpoint)
    // {
    //     vtkgpoint = dataset->AllocatePointGhostArray();
    // }

    UList<uint8_t> gcell = vtk::Tools::asUList(vtkgcell, cellMap.size());

    vtkgcell->FillValue(0); // Initialize to zero

    forAll(cellMap, i)
    {
        const label cellType = types[cellMap[i]];

        if (cellType == cellCellStencil::HOLE)
        {
            // Need duplicate (not just HIDDENCELL) for it to be properly
            // recognized
            gcell[i] |=
            (
                vtkDataSetAttributes::DUPLICATECELL
              | vtkDataSetAttributes::HIDDENCELL
            );
        }
        #if 0
        // No special treatment for INTERPOLATED.
        // This causes duplicate/overlapping values, but avoids holes
        // in the results
        else if (cellType == cellCellStencil::INTERPOLATED)
        {
            gcell[i] |=
            (
                vtkDataSetAttributes::DUPLICATECELL
            );
        }
        #endif
    }

    dataset->GetCellData()->AddArray(vtkgcell);
}


void Foam::vtk::fvMeshAdaptor::applyGhostingBoundary(const labelUList& types)
{
    // BOUNDARY

    if (types.empty() || !usingBoundary())
    {
        return;
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const label patchId : patchIds_)
    {
        const polyPatch& pp = patches[patchId];
        const word& longName = pp.name();

        auto iter = cachedVtp_.find(longName);
        if (!iter.good() || !iter.val().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            Pout<<"Cache miss for VTP patch " << longName << endl;
            continue;
        }

        foamVtpData& vtpData = iter.val();
        auto dataset = vtpData.dataset;

        auto vtkgcell = dataset->GetCellGhostArray();
        if (!vtkgcell)
        {
            vtkgcell = dataset->AllocateCellGhostArray();
        }


        // Determine face ghosting based on interior cells
        const labelUList& bCells = pp.faceCells();

        const label len = bCells.size();

        UList<uint8_t> gcell = vtk::Tools::asUList(vtkgcell, len);

        vtkgcell->FillValue(0); // Initialize to zero

        for (label i=0; i < len; ++i)
        {
            const label celli = bCells[i];

            const label cellType = types[celli];

            if (cellType == cellCellStencil::HOLE)
            {
                // Need duplicate (not just HIDDENCELL) for it to be properly
                // recognized
                gcell[i] |=
                (
                    vtkDataSetAttributes::DUPLICATECELL
                  | vtkDataSetAttributes::HIDDENCELL
                );
            }
            #if 0
            // No special treatment for INTERPOLATED.
            // This causes duplicate/overlapping values, but avoids holes
            // in the results
            else if (cellType == cellCellStencil::INTERPOLATED)
            {
                gcell[i] |=
                (
                    vtkDataSetAttributes::DUPLICATECELL
                );
            }
            #endif
        }

        dataset->GetCellData()->AddArray(vtkgcell);
    }
}


// These need to be rebuild here, since the component mesh may just have point
// motion without topology changes.
void Foam::vtk::fvMeshAdaptor::applyGhosting()
{
    const auto* stencilPtr =
        mesh_.cfindObject<cellCellStencilObject>
        (
            cellCellStencilObject::typeName
        );

    if (stencilPtr)
    {
        const labelUList& types = stencilPtr->cellTypes();

        applyGhostingInternal(types);
        applyGhostingBoundary(types);
    }
}


// ************************************************************************* //
