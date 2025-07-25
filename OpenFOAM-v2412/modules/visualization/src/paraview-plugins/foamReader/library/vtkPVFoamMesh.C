/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "faMesh.H"
#include "fvMeshSubset.H"
#include "vtkPVFoamReader.h"
#include "uindirectPrimitivePatch.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume()
{
    // Convert internalMesh - actually only one part
    const List<label> partIds =
        rangeVolume_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        foamVtuData& vtuData = cachedVtu_(longName);

        vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
        if (vtuData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo << "reuse " << longName << nl;
                vtuData.reuse(); // reuse
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                DebugInfo << "move points " << longName << nl;
                vtkgeom = vtuData.getCopy();
                vtkgeom->SetPoints(vtuData.points(mesh));
            }
        }

        if (!vtkgeom)
        {
            // Nothing usable from cache - create new geometry
            vtkgeom = vtuData.internal(mesh, this->decomposePoly_);
        }

        vtuData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshLagrangian()
{
    const List<label> partIds =
        rangeClouds_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word cloudName = getFoamName(longName);

        // Direct conversion, no caching for Lagrangian
        cachedVtp_(longName).set
        (
            lagrangianVTKMesh(mesh, cloudName)
        );
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPatches()
{
    const List<label> partIds =
        rangePatches_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word partName = getFoamName(longName);
        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                DebugInfo << "reuse " << longName << nl;
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Point movement on single patch is OK

                const labelUList& patchIds = vtpData.additionalIds();
                if (patchIds.size() == 1)
                {
                    vtkgeom = vtpData.getCopy();
                    vtkgeom->SetPoints
                    (
                        vtk::Tools::Patch::points(patches[patchIds[0]])
                    );
                    continue;
                }
            }
        }

        if (longName.starts_with("group/"))
        {
            // Patch group. Collect patch faces.

            vtpData.clear(); // Remove any old mappings

            const labelList& patchIds =
                patches.groupPatchIDs().lookup(partName, labelList());

            DebugInfo
                << "Creating VTK mesh for patches [" << patchIds <<"] "
                << longName << nl;

            // Store good patch ids as additionalIds
            vtpData.additionalIds().reserve(patchIds.size());

            label sz = 0;
            for (const auto id : patchIds)
            {
                const label n = patches[id].size();

                if (n)
                {
                    vtpData.additionalIds().append(id);
                    sz += n;
                }
            }
            Foam::sort(vtpData.additionalIds());

            // Temporarily (mis)use cellMap for face labels
            DynamicList<label>& faceLabels = vtpData.cellMap();
            faceLabels.reserve(sz);

            for (const auto id : vtpData.additionalIds())
            {
                const auto& pp = patches[id];
                forAll(pp, i)
                {
                    faceLabels.append(pp.start()+i);
                }
            }

            if (faceLabels.size())
            {
                uindirectPrimitivePatch pp
                (
                    UIndirectList<face>(mesh.faces(), faceLabels),
                    mesh.points()
                );

                vtkgeom = vtk::Tools::Patch::mesh(pp);
            }

            faceLabels.clear();  // Unneeded
        }
        else
        {
            vtpData.clear(); // Remove any old mappings

            const label patchId = patches.findPatchID(partName);

            DebugInfo
                << "Creating VTK mesh for patch [" << patchId <<"] "
                << partName << nl;

            if (patchId >= 0)
            {
                // Store good patch id as additionalIds
                vtpData.additionalIds() = {patchId};

                vtkgeom = vtk::Tools::Patch::mesh(patches[patchId]);
            }
        }

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            cachedVtp_.erase(longName);
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshArea()
{
    // Convert areaMesh - actually only one part
    const List<label> partIds =
        rangeArea_.intersection(selectedPartIds_);

    if (partIds.empty())
    {
        return;
    }
    else if (!areaMeshPtr_ && volMeshPtr_)
    {
        areaMeshPtr_ = new faMesh(*volMeshPtr_);
    }

    if (!areaMeshPtr_)
    {
        return;
    }

    const faMesh& mesh = *areaMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;

// This needs checking
// #if 0
//         if (vtpData.nPoints())
//         {
//             if (meshState_ == polyMesh::UNCHANGED)
//             {
//                 DebugInfo << "reuse " << longName << nl;
//                 vtpData.reuse(); // reuse
//                 continue;
//             }
//             else if (meshState_ == polyMesh::POINTS_MOVED)
//             {
//                 DebugInfo << "move points " << longName << nl;
//                 vtkgeom = vtpData.getCopy();
//                 vtkgeom->SetPoints(movePoints(mesh, vtpData));
//             }
//         }
// #endif

        if (!vtkgeom)
        {
            vtpData.clear(); // Remove any old mappings

            // Nothing usable from cache - create new geometry
            vtkgeom = vtk::Tools::Patch::mesh(mesh.patch());
        }

        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellZones()
{
    const List<label> partIds =
        rangeCellZones_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (partIds.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    const cellZoneMesh& zMesh = mesh.cellZones();
    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word zoneName = getFoamName(longName);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }

        DebugInfo
            << "Creating VTK mesh for cellZone[" << zoneId << "] "
            << zoneName << nl;

        foamVtuData& vtuData = cachedVtu_(longName);

        vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
        if (vtuData.nPoints() && vtuData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo << "reuse " << longName << nl;
                vtuData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                DebugInfo << "move points " << longName << nl;
                vtkgeom = vtuData.getCopy();
                vtkgeom->SetPoints(vtuData.points(mesh, vtuData.pointMap()));
            }
        }

        if (!vtkgeom)
        {
            fvMeshSubset subsetter(mesh, zMesh[zoneId]);

            vtkgeom = vtuData.subset(subsetter, this->decomposePoly_);
        }

        vtuData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellSets()
{
    const List<label> partIds =
        rangeCellSets_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word partName = getFoamName(longName);

        DebugInfo << "Creating VTK mesh for cellSet=" << partName << nl;

        foamVtuData& vtuData = cachedVtu_(longName);

        vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
        if (vtuData.nPoints() && vtuData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo << "reuse " << longName << nl;
                vtuData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                DebugInfo << "move points " << longName << nl;
                vtkgeom = vtuData.getCopy();
                vtkgeom->SetPoints(vtuData.points(mesh, vtuData.pointMap()));
            }
        }

        if (!vtkgeom)
        {
            fvMeshSubset subsetter(mesh, cellSet(mesh, partName));

            vtkgeom = vtuData.subset(subsetter, this->decomposePoly_);
        }

        vtuData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones()
{
    const List<label> partIds =
        rangeFaceZones_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (partIds.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    const faceZoneMesh& zMesh = mesh.faceZones();
    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word zoneName = getFoamName(longName);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }
        DebugInfo
            << "Creating VTKmesh for faceZone[" << zoneId << "] "
            << zoneName << nl;

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                DebugInfo <<"reuse " << longName << nl;
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Need point maps etc - not worth it at the moment
            }
        }

        if (!vtkgeom)
        {
            vtpData.clear();  // No additional ids, maps

            const primitiveFacePatch& pp = zMesh[zoneId]();
            vtkgeom = vtk::Tools::Patch::mesh(pp);
        }

        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets()
{
    const List<label> partIds =
        rangeFaceSets_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word partName = getFoamName(longName);

        DebugInfo << "Creating VTK mesh for faceSet=" << partName << nl;

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                DebugInfo <<"reuse " << longName << nl;
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Need point maps etc - not worth it at the moment
            }
        }

        if (!vtkgeom)
        {
            vtpData.clear();  // No other additional ids, maps

            // Misuse cellMap for face labels - sorted order for reliability
            vtpData.cellMap() = faceSet(mesh, partName).sortedToc();

            if (vtpData.cellMap().size())
            {
                uindirectPrimitivePatch pp
                (
                    UIndirectList<face>(mesh.faces(), vtpData.cellMap()),
                    mesh.points()
                );

                vtkgeom = vtk::Tools::Patch::mesh(pp);
            }
        }

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            cachedVtp_.erase(longName);
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones()
{
    const List<label> partIds =
        rangePointZones_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (partIds.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    const pointZoneMesh& zMesh = mesh.pointZones();
    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word zoneName = getFoamName(longName);
        const label zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints() && vtpData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo << "reusing " << longName << nl;
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                DebugInfo << "move points " << longName << nl;
                vtkgeom = vtpData.getCopy();
            }
        }

        if (!vtkgeom)
        {
            // First time, or topo change
            vtkgeom = vtkSmartPointer<vtkPolyData>::New();
            vtpData.pointMap() = zMesh[zoneId];
        }

        const pointField& points = mesh.points();
        const labelUList& pointMap = vtpData.pointMap();

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

        vtkpoints->SetNumberOfPoints(pointMap.size());
        forAll(pointMap, pointi)
        {
            vtkpoints->SetPoint(pointi, points[pointMap[pointi]].v_);
        }

        vtkgeom->SetPoints(vtkpoints);
        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointSets()
{
    const List<label> partIds =
        rangePointSets_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word partName = getFoamName(longName);

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints() && vtpData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                DebugInfo << "reusing " << longName << nl;
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                DebugInfo << "move points " << longName << nl;
                vtkgeom = vtpData.getCopy();
            }
        }

        if (!vtkgeom)
        {
            // First time, or topo change
            vtkgeom = vtkSmartPointer<vtkPolyData>::New();
            vtpData.pointMap() = pointSet(mesh, partName).sortedToc();
        }

        const pointField& points = mesh.points();
        const labelUList& pointMap = vtpData.pointMap();

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

        vtkpoints->SetNumberOfPoints(pointMap.size());
        forAll(pointMap, pointi)
        {
            vtkpoints->SetPoint(pointi, points[pointMap[pointi]].v_);
        }

        vtkgeom->SetPoints(vtkpoints);
        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


// ************************************************************************* //
