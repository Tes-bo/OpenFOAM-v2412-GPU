/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "foamVtkCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::fileExtension
({
    { fileTag::POLY_DATA, "vtp" },
    { fileTag::UNSTRUCTURED_GRID, "vtu" },
    { fileTag::MULTI_BLOCK, "vtm" },
    // { fileTag::COLLECTION, "pvd" },
});


const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::fileContentVersions
({
    { fileTag::POLY_DATA, "0.1" },
    { fileTag::UNSTRUCTURED_GRID, "0.1" },
    { fileTag::MULTI_BLOCK, "1.0" },
    // { fileTag::COLLECTION, "0.1" },
});


const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::fileTagNames
({
    { fileTag::VTK_FILE, "VTKFile" },
    { fileTag::DATA_ARRAY, "DataArray" },
    { fileTag::BLOCK, "Block" },
    { fileTag::PIECE, "Piece" },
    { fileTag::DATA_SET, "DataSet" },
    { fileTag::POINTS, "Points" },
    { fileTag::CELLS, "Cells" },
    { fileTag::POLYS, "Polys" },
    { fileTag::VERTS, "Verts" },
    { fileTag::LINES, "Lines" },
    { fileTag::CELL_DATA, "CellData" },
    { fileTag::POINT_DATA, "PointData" },
    { fileTag::FIELD_DATA, "FieldData" },
    { fileTag::POLY_DATA, "PolyData" },
    { fileTag::UNSTRUCTURED_GRID, "UnstructuredGrid" },
    { fileTag::MULTI_BLOCK, "vtkMultiBlockDataSet" },
    // { fileTag::COLLECTION, "Collection" },
});


const Foam::Enum
<
    Foam::vtk::fileAttr
>
Foam::vtk::fileAttrNames
({
    { fileAttr::OFFSET, "offset" },
    { fileAttr::NUMBER_OF_COMPONENTS, "NumberOfComponents" },
    { fileAttr::NUMBER_OF_TUPLES, "NumberOfTuples" },
    { fileAttr::NUMBER_OF_POINTS, "NumberOfPoints" },
    { fileAttr::NUMBER_OF_CELLS, "NumberOfCells" },
    { fileAttr::NUMBER_OF_POLYS, "NumberOfPolys" },
    { fileAttr::NUMBER_OF_VERTS, "NumberOfVerts" },
    { fileAttr::NUMBER_OF_LINES, "NumberOfLines" },
});


const Foam::Enum
<
    Foam::vtk::dataArrayAttr
>
Foam::vtk::dataArrayAttrNames
({
    { dataArrayAttr::POINTS, "Points" },
    { dataArrayAttr::OFFSETS, "offsets" },
    { dataArrayAttr::CONNECTIVITY, "connectivity" },
    { dataArrayAttr::TYPES, "types" },
    { dataArrayAttr::FACES, "faces" },
    { dataArrayAttr::FACEOFFSETS, "faceoffsets" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Legacy

const Foam::word Foam::vtk::legacy::fileExtension("vtk");

const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::legacy::contentNames
({
    { vtk::fileTag::POLY_DATA, "POLYDATA" },
    { vtk::fileTag::UNSTRUCTURED_GRID, "UNSTRUCTURED_GRID" },
});


const Foam::Enum
<
    Foam::vtk::fileTag
>
Foam::vtk::legacy::fileTagNames
({
    { vtk::fileTag::POINTS, "POINTS" },
    { vtk::fileTag::CELLS, "CELLS" },
    { vtk::fileTag::POLYS, "POLYGONS" },
    { vtk::fileTag::VERTS, "VERTICES" },
    { vtk::fileTag::LINES, "LINES" },
    { vtk::fileTag::CELL_DATA, "CELL_DATA" },
    { vtk::fileTag::POINT_DATA, "POINT_DATA" },
});


const Foam::Enum
<
    Foam::vtk::dataArrayAttr
>
Foam::vtk::legacy::dataArrayAttrNames
({
    { vtk::dataArrayAttr::OFFSETS, "OFFSETS" },
    { vtk::dataArrayAttr::CONNECTIVITY, "CONNECTIVITY" },
});


// ************************************************************************* //
