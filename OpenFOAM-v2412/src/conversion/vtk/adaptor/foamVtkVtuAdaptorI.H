/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pts.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const point& p : pts)
    {
        vtkpoints->SetPoint(pointId++, p.v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh,
    const labelUList& pointMap
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pointMap.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const label meshPointi : pointMap)
    {
        vtkpoints->SetPoint(pointId++, pts[meshPointi].v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkUnstructuredGrid>
Foam::vtk::vtuAdaptor::internal
(
    const fvMesh& mesh,
    const bool decompPoly
)
{
    const vtk::vtuSizing::contentType output
    (
        #ifdef VTK_CELL_ARRAY_V2
        vtk::vtuSizing::contentType::INTERNAL2
        #else
        vtk::vtuSizing::contentType::INTERNAL1
        #endif
    );

    vtk::vtuSizing sizing(mesh, decompPoly);

    auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();

    UList<uint8_t> cellTypesUL
    (
        vtk::Tools::asUList(cellTypes, sizing.nFieldCells())
    );

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkIdTypeArray>::New();

    auto cellOffsets   = vtkSmartPointer<vtkIdTypeArray>::New();
    auto faceLocations = vtkSmartPointer<vtkIdTypeArray>::New();

    const auto nConnect
    (
        sizing.sizeOf(output, vtk::vtuSizing::slotType::CELLS)
    );

    UList<vtkIdType> cellOffsetsUL
    (
        vtk::Tools::asUList
        (
            cellOffsets,
            sizing.sizeOf(output, vtk::vtuSizing::slotType::CELLS_OFFSETS)
        )
    );

    #ifdef VTK_CELL_ARRAY_V2

    auto cellConnect = vtkSmartPointer<vtkIdTypeArray>::New();

    UList<vtkIdType> cellsUL
    (
        vtk::Tools::asUList(cellConnect, nConnect)
    );

    #else

    cells->GetData()->SetNumberOfTuples(sizing.nFieldCells());

    UList<vtkIdType> cellsUL
    (
        cells->WritePointer(sizing.nFieldCells(), nConnect),
        nConnect
    );

    #endif

    UList<vtkIdType> facesUL
    (
        vtk::Tools::asUList
        (
            faces,
            sizing.sizeOf(output, vtk::vtuSizing::slotType::FACES)
        )
    );

    UList<vtkIdType> faceLocationsUL
    (
        vtk::Tools::asUList
        (
            faceLocations,
            sizing.sizeOf(output, vtk::vtuSizing::slotType::FACES_OFFSETS)
        )
    );


    sizing.populateInternal
    (
        mesh,
        cellTypesUL,
        cellsUL, cellOffsetsUL,
        facesUL, faceLocationsUL,
        static_cast<foamVtkMeshMaps&>(*this),
        output
    );


    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkmesh->SetPoints(this->points(mesh));

    #ifdef VTK_CELL_ARRAY_V2

    // Move into a vtkCellArray
    cells->SetData(cellOffsets, cellConnect);

    if (facesUL.size())
    {
        vtkmesh->SetCells(cellTypes, cells, faceLocations, faces);
    }
    else
    {
        vtkmesh->SetCells(cellTypes, cells, nullptr, nullptr);
    }
    #else

    if (facesUL.size())
    {
        vtkmesh->SetCells(cellTypes, cellOffsets, cells, faceLocations, faces);
    }
    else
    {
        vtkmesh->SetCells(cellTypes, cellOffsets, cells, nullptr, nullptr);
    }

    #endif

    return vtkmesh;
}


// ************************************************************************* //
