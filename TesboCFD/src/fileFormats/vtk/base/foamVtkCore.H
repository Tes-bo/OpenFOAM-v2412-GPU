/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Namespace
    Foam::vtk

Description
    Namespace for handling VTK input/output.

SourceFiles
    foamVtkCore.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_core_H
#define Foam_vtk_core_H

#include <cstdint>
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Enumerations

    //- The context when outputting a VTK file (XML or legacy).
    enum OutputContext : uint8_t
    {
        INLINE,     //!< Generate header and inline data
        HEADER,     //!< Generate header only
        APPEND      //!< Generate append-data
    };


    //- The output format type for file contents.
    //  Upper bits for output type, lower bits for the format itself.
    enum class formatType : uint8_t
    {
        INLINE_ASCII  = 0,    //!< XML inline ASCII, asciiFormatter
        INLINE_BASE64 = 0x01, //!< XML inline base64, base64Formatter
        APPEND_BASE64 = 0x11, //!< XML append base64, appendBase64Formatter
        APPEND_BINARY = 0x12, //!< XML append raw binary, appendRawFormatter
        LEGACY_ASCII  = 0x20, //!< Legacy ASCII, legacyAsciiFormatter
        LEGACY_BINARY = 0x22, //!< Legacy raw binary, legacyRawFormatter
    };

    //- Test for vtk append format (xml)
    inline bool isAppend(enum formatType fmt) noexcept
    {
        return (uint8_t(fmt) & 0x10);
    }

    //- Test for vtk legacy format
    inline bool isLegacy(enum formatType fmt) noexcept
    {
        return (uint8_t(fmt) & 0x20);
    }


    //- Equivalent to enumeration in "vtkCellType.h" (should be uint8_t)
    enum cellType : uint8_t
    {
        VTK_EMPTY_CELL       = 0,
        VTK_VERTEX           = 1,
        VTK_POLY_VERTEX      = 2,
        VTK_LINE             = 3,
        VTK_POLY_LINE        = 4,
        VTK_TRIANGLE         = 5,
        VTK_TRIANGLE_STRIP   = 6,
        VTK_POLYGON          = 7,
        VTK_PIXEL            = 8,
        VTK_QUAD             = 9,
        VTK_TETRA            = 10,
        VTK_VOXEL            = 11,
        VTK_HEXAHEDRON       = 12,
        VTK_WEDGE            = 13,
        VTK_PYRAMID          = 14,
        VTK_PENTAGONAL_PRISM = 15,
        VTK_HEXAGONAL_PRISM  = 16,
        VTK_POLYHEDRON       = 42
    };


    //- Some common XML tags for vtk files
    enum class fileTag : uint8_t
    {
        VTK_FILE,               //!< "VTKFile"
        DATA_ARRAY,             //!< "DataArray"
        BLOCK,                  //!< "Block"
        PIECE,                  //!< "Piece"
        DATA_SET,               //!< "DataSet"
        POINTS,                 //!< "Points"
        CELLS,                  //!< "Cells"
        POLYS,                  //!< "Polys"
        VERTS,                  //!< "Verts"
        LINES,                  //!< "Lines"
        CELL_DATA,              //!< "CellData"
        POINT_DATA,             //!< "PointData"
        FIELD_DATA,             //!< "FieldData"
        POLY_DATA,              //!< "PolyData"
        UNSTRUCTURED_GRID,      //!< "UnstructuredGrid"
        MULTI_BLOCK,            //!< "vtkMultiBlockDataSet"
    };

    //- File extension (without ".") for some vtk XML file content types
    extern const Foam::Enum<fileTag> fileExtension;

    //- Version string for some vtk XML file content types
    extern const Foam::Enum<fileTag> fileContentVersions;

    //- Strings corresponding to the vtk XML tags
    extern const Foam::Enum<fileTag> fileTagNames;

    //- Some common XML attributes for vtk files
    enum class fileAttr : uint8_t
    {
        OFFSET,                 //!< "offset"
        NUMBER_OF_COMPONENTS,   //!< "NumberOfComponents"
        NUMBER_OF_TUPLES,       //!< "NumberOfTuples"
        NUMBER_OF_POINTS,       //!< "NumberOfPoints"
        NUMBER_OF_CELLS,        //!< "NumberOfCells"
        NUMBER_OF_POLYS,        //!< "NumberOfPolys"
        NUMBER_OF_VERTS,        //!< "NumberOfVerts"
        NUMBER_OF_LINES,        //!< "NumberOfLines"
    };

    //- Strings corresponding to the vtk XML attributes
    extern const Foam::Enum<fileAttr> fileAttrNames;

    //- Some common names for XML DataArray entries
    enum class dataArrayAttr : uint8_t
    {
        POINTS,                 //!< "Points"
        OFFSETS,                //!< "offsets"
        CONNECTIVITY,           //!< "connectivity"
        TYPES,                  //!< "types"
        FACES,                  //!< "faces"
        FACEOFFSETS,            //!< "faceoffsets"
    };

    //- Strings corresponding to the vtk XML DataArray attributes
    extern const Foam::Enum<dataArrayAttr> dataArrayAttrNames;


/*---------------------------------------------------------------------------*\
                               Namespace legacy
\*---------------------------------------------------------------------------*/

namespace legacy
{

    //- Legacy file extension ("vtk")
    extern const word fileExtension;

    //- Legacy content names (POLYDATA, UNSTRUCTURED_GRID)
    extern const Foam::Enum<vtk::fileTag> contentNames;

    //- Legacy file tags (eg, LINES, CELL_DATA, POINT_DATA, ...)
    extern const Foam::Enum<vtk::fileTag> fileTagNames;

    //- Legacy attributes (eg, OFFSETS)
    extern const Foam::Enum<dataArrayAttr> dataArrayAttrNames;


} // End namespace legacy


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
