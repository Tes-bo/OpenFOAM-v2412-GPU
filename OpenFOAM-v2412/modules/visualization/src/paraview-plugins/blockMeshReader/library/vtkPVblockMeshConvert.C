/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "vtkPVblockMesh.H"
#include "vtkPVblockMeshReader.h"

// OpenFOAM includes
#include "blockMesh.H"
#include "Time.H"
#include "foamVtkTools.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVblockMesh::convertMeshBlocks
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    DebugInfo << "<beg> convertMeshBlocks" << endl;

    const Map<string> blockStatus = getSelectedArrayMap
    (
        reader_->GetBlockSelection()
    );

    arrayRange& range = rangeBlocks_;
    range.block(blockNo);   // set output block
    label datasetNo = 0;    // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const pointField blkPoints(blkMesh.vertices(true));

    vtkIdType nodeIds[8];  // Exact space for VTK_HEXAHEDRON vertices
    int blockId = -1;
    for (auto partId : range)
    {
        ++blockId; // Increment first
        if (!blockStatus.contains(partId))
        {
            continue;
        }
        const auto& longName = blockStatus[partId];

        const blockDescriptor& blockDef = blkMesh[blockId];
        const labelList& blockLabels = blockDef.blockShape();

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkpoints->SetNumberOfPoints(blockLabels.size());

        forAll(blockLabels, pointi)
        {
            vtkpoints->SetPoint
            (
                pointi,
                blkPoints[blockLabels[pointi]].v_
            );
            nodeIds[pointi] = pointi;
        }

        auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkmesh->Allocate(1);
        vtkmesh->InsertNextCell
        (
            VTK_HEXAHEDRON,
            8,
            nodeIds
        );

        vtkmesh->SetPoints(vtkpoints);

        addToBlock(output, vtkmesh, range, datasetNo, longName);
        ++datasetNo;
    }


    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    DebugInfo << "<end> convertMeshBlocks" << endl;
}


void Foam::vtkPVblockMesh::convertMeshEdges
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    const Map<string> edgeStatus = getSelectedArrayMap
    (
        reader_->GetCurvedEdgesSelection()
    );

    arrayRange& range = rangeEdges_;

    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const blockEdgeList& edges = blkMesh.edges();

    // Use 8 uniform divisions for visualization purposes (smooth enough)
    constexpr label nEdgeDivs = 4;
    List<vtkIdType> pointIds(nEdgeDivs + 1);

    // Edge point and weighting factors for a specified hex edge
    pointField cedgePoints;
    scalarList cedgeWeights;

    int edgeId = -1;
    for (auto partId : range)
    {
        ++edgeId; // Increment first
        if (!edgeStatus.contains(partId))
        {
            continue;
        }
        const auto& longName = edgeStatus[partId];

        // Search each block
        for (const blockDescriptor& blockDef : blkMesh)
        {
            edgeList blkEdges(blockDef.blockShape().edges());

            // Find the corresponding edge within the block
            label foundEdgei = -1;
            forAll(blkEdges, blkEdgei)
            {
                if (edges[edgeId].compare(blkEdges[blkEdgei]))
                {
                    foundEdgei = blkEdgei;
                    break;
                }
            }

            if (foundEdgei != -1)
            {
                blockDef.edgePointsWeights
                (
                    foundEdgei,
                    cedgePoints,
                    cedgeWeights,
                    nEdgeDivs
                );

                const pointField edgePoints
                (
                    blkMesh.globalPosition(cedgePoints)
                );

                auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
                vtkpoints->SetNumberOfPoints(edgePoints.size());

                pointIds.resize(edgePoints.size());

                forAll(edgePoints, pointi)
                {
                    vtkpoints->SetPoint(pointi, edgePoints[pointi].v_);
                    pointIds[pointi] = pointi;
                }

                auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

                vtkmesh->Allocate(1);
                vtkmesh->InsertNextCell
                (
                    VTK_POLY_LINE,
                    edgePoints.size(),
                    pointIds.cdata()
                );

                vtkmesh->SetPoints(vtkpoints);

                addToBlock(output, vtkmesh, range, datasetNo, longName);
                ++datasetNo;

                break;
            }
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    DebugInfo << "<end> convertMeshEdges" << endl;
}


void Foam::vtkPVblockMesh::convertMeshCorners
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = rangeCorners_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const pointField blkPoints(blkMesh.vertices(true));

    DebugInfo << "<beg> " << FUNCTION_NAME << endl;

    if (true)  // Or some flag or other condition
    {
        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkpoints->SetNumberOfPoints(blkPoints.size());

        forAll(blkPoints, pointi)
        {
            vtkpoints->SetPoint(pointi, blkPoints[pointi].v_);
        }

        auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

        vtkmesh->SetPoints(vtkpoints);
        vtkmesh->SetVerts(vtk::Tools::identityVertices(blkPoints.size()));

        addToBlock(output, vtkmesh, range, datasetNo, range.name());
        ++datasetNo;
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    DebugInfo << "<end> " << FUNCTION_NAME << endl;
}


// ************************************************************************* //
