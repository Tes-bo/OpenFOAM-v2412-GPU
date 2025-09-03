/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModel.H"
#include "vectorIOField.H"
#include "triangle.H"
#include "stringOps.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(vtkUnstructuredReader, 1);
}

const Foam::Enum
<
    Foam::vtkUnstructuredReader::vtkDataType
>
Foam::vtkUnstructuredReader::vtkDataTypeNames
({
    { vtkDataType::VTK_INT, "int" },
    // Not yet required: { vtkDataType::VTK_INT64, "vtktypeint64" },
    { vtkDataType::VTK_UINT, "unsigned_int" },
    { vtkDataType::VTK_LONG, "long" },
    { vtkDataType::VTK_ULONG, "unsigned_long" },
    { vtkDataType::VTK_FLOAT, "float" },
    { vtkDataType::VTK_DOUBLE, "double" },
    { vtkDataType::VTK_STRING, "string" },
    { vtkDataType::VTK_ID, "vtkIdType" },
});


const Foam::Enum
<
    Foam::vtkUnstructuredReader::vtkDataSetType
>
Foam::vtkUnstructuredReader::vtkDataSetTypeNames
({
    { vtkDataSetType::VTK_FIELD, "FIELD" },
    { vtkDataSetType::VTK_SCALARS, "SCALARS" },
    { vtkDataSetType::VTK_VECTORS, "VECTORS" },
});


const Foam::Enum
<
    Foam::vtkUnstructuredReader::parseMode
>
Foam::vtkUnstructuredReader::parseModeNames
({
    { parseMode::NOMODE, "NOMODE" },
    { parseMode::UNSTRUCTURED_GRID, "UNSTRUCTURED_GRID" },
    { parseMode::POLYDATA, "POLYDATA" },
    { parseMode::CELL_DATA, "CELL_DATA" },
    { parseMode::POINT_DATA, "POINT_DATA" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read N elements into List
template<class T>
static inline void readBlock(Istream& is, const label n, List<T>& list)
{
    list.resize(n);
    for (T& val : list)
    {
        is >> val;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::vtkUnstructuredReader::warnUnhandledType
(
    const Istream& is,
    const label type,
    labelHashSet& warningGiven
)
{
    if (warningGiven.insert(type))
    {
        IOWarningInFunction(is)
            << "Skipping unknown cell type " << type << nl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkUnstructuredReader::readOffsetsConnectivity
(
    ISstream& is,
    const char* entryName,
    const label nOffsets,
    labelList& offsets,
    const label nConnectivity,
    labelList& connectivity
)
{
    token tok;

    is.read(tok);
    if (!tok.isWord("OFFSETS"))
    {
        FatalIOErrorInFunction(is)
            << "Expected OFFSETS for " << entryName
            << ", found "
            << tok.info() << nl
            << exit(FatalIOError);
    }
    is.getLine(nullptr);  // Consume rest of line
    readBlock(is, nOffsets, offsets);

    is.read(tok);
    if (!tok.isWord("CONNECTIVITY"))
    {
        FatalIOErrorInFunction(is)
            << "Expected CONNECTIVITY for " << entryName
            << ", found " << tok.info() << nl
            << exit(FatalIOError);
    }
    is.getLine(nullptr);  // Consume rest of line
    readBlock(is, nConnectivity, connectivity);
}


void Foam::vtkUnstructuredReader::extractCells
(
    const Istream& is,
    const labelUList& cellTypes,
    const labelUList& elemOffsets,
    const labelUList& elemVerts
)
{
    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& pyr = cellModel::ref(cellModel::PYR);
    const cellModel& tet = cellModel::ref(cellModel::TET);

    // Pass 0: sizing
    label nCells = 0, nFaces = 0, nLines = 0;
    forAll(cellTypes, elemi)
    {
        switch (cellTypes[elemi])
        {
            case vtk::cellType::VTK_LINE:
            case vtk::cellType::VTK_POLY_LINE:
            {
                ++nLines;
            }
            break;

            case vtk::cellType::VTK_TRIANGLE:
            case vtk::cellType::VTK_QUAD:
            case vtk::cellType::VTK_POLYGON:
            {
                ++nFaces;
            }
            break;

            case vtk::cellType::VTK_TETRA:
            case vtk::cellType::VTK_PYRAMID:
            case vtk::cellType::VTK_WEDGE:
            case vtk::cellType::VTK_HEXAHEDRON:
            {
                ++nCells;
            }
            break;

            default:
                break;
        }
    }

    label celli = cells_.size();
    cells_.resize(celli + nCells);
    cellMap_.resize(cells_.size(), -1);

    label facei = faces_.size();
    faces_.resize(facei + nFaces);
    faceMap_.resize(faces_.size(), -1);

    label linei = lines_.size();
    lines_.resize(linei + nLines);
    lineMap_.resize(lines_.size(), -1);

    // General scratch space for cell vertices
    labelList cellPoints(16);

    label dataIndex = 0;  // Addressing into vertices stream

    // To mark whether unhandled type has been visited.
    labelHashSet warningGiven;


    forAll(cellTypes, elemi)
    {
        // Vertices per element - from offsets or embedded size
        const label nVerts =
        (
            elemOffsets.empty()
          ? elemVerts[dataIndex++]
          : (elemOffsets[elemi+1] - elemOffsets[elemi])
        );

        // The addresseed vertices
        const labelSubList verts(elemVerts, nVerts, dataIndex);
        dataIndex += nVerts;

        switch (cellTypes[elemi])
        {
            case vtk::cellType::VTK_VERTEX:
            {
                warnUnhandledType(is, cellTypes[elemi], warningGiven);
                if (nVerts != 1)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 1 for VTK_VERTEX, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
            }
            break;

            case vtk::cellType::VTK_POLY_VERTEX:
            {
                warnUnhandledType(is, cellTypes[elemi], warningGiven);
            }
            break;

            case vtk::cellType::VTK_LINE:
            {
                if (nVerts != 2)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 2 for VTK_LINE, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                lineMap_[linei] = elemi;

                // Same vertex ordering
                lines_[linei++] = verts;
            }
            break;

            case vtk::cellType::VTK_POLY_LINE:
            {
                lineMap_[linei] = elemi;

                // Same vertex ordering
                lines_[linei++] = verts;
            }
            break;

            case vtk::cellType::VTK_TRIANGLE:
            {
                if (nVerts != 3)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 3 for VTK_TRIANGLE, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                faceMap_[facei] = elemi;

                // Same vertex ordering
                static_cast<labelList&>(faces_[facei++]) = verts;
            }
            break;

            case vtk::cellType::VTK_QUAD:
            {
                if (nVerts != 4)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 4 for VTK_QUAD, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                faceMap_[facei] = elemi;

                // Same vertex ordering
                static_cast<labelList&>(faces_[facei++]) = verts;
            }
            break;

            case vtk::cellType::VTK_POLYGON:
            {
                faceMap_[facei] = elemi;

                // Same vertex ordering
                static_cast<labelList&>(faces_[facei++]) = verts;
            }
            break;

            case vtk::cellType::VTK_TETRA:
            {
                if (nVerts != 4)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 4 for VTK_TETRA, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                cellMap_[celli] = elemi;

                // Same vertex ordering
                cells_[celli++].reset(tet, verts, true);
            }
            break;

            case vtk::cellType::VTK_PYRAMID:
            {
                if (nVerts != 5)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 5 for VTK_PYRAMID, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                cellMap_[celli] = elemi;

                // Same vertex ordering
                cells_[celli++].reset(pyr, verts, true);
            }
            break;

            case vtk::cellType::VTK_WEDGE:
            {
                if (nVerts != 6)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 6 for VTK_WEDGE, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                cellMap_[celli] = elemi;

                // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
                labelSubList shape(cellPoints, nVerts);
                shape[0] = verts[0];
                shape[2] = verts[1];
                shape[1] = verts[2];
                shape[3] = verts[3];
                shape[5] = verts[4];
                shape[4] = verts[5];

                cells_[celli++].reset(prism, shape, true);
            }
            break;

            case vtk::cellType::VTK_HEXAHEDRON:
            {
                if (nVerts != 8)
                {
                    FatalIOErrorInFunction(is)
                        << "Expected size 8 for VTK_HEXAHEDRON, found "
                        << nVerts << nl
                        << exit(FatalIOError);
                }
                cellMap_[celli] = elemi;

                // Same vertex ordering
                cells_[celli++].reset(hex, verts, true);
            }
            break;

            default:
            {
                warnUnhandledType(is, cellTypes[elemi], warningGiven);
            }
            break;
        }
    }

    DebugInfo
        << "Read"
        << " cells:" << celli
        << " faces:" << facei
        << " lines:" << linei
        << nl;

    cells_.resize(celli);
    cellMap_.resize(celli);

    faces_.resize(facei);
    faceMap_.resize(facei);

    lines_.resize(linei);
    lineMap_.resize(linei);
}


void Foam::vtkUnstructuredReader::readField
(
    ISstream& inFile,
    objectRegistry& obj,
    const word& arrayName,
    const word& dataType,
    const label size
) const
{
    if (vtkDataTypeNames.found(dataType))
    {
        switch (vtkDataTypeNames[dataType])
        {
            case VTK_INT:
            case VTK_INT64:
            case VTK_UINT:
            case VTK_LONG:
            case VTK_ULONG:
            case VTK_ID:
            {
                auto fieldVals = autoPtr<labelIOField>::New
                (
                    IOobject(arrayName, "", obj),
                    size
                );
                readBlock(inFile, fieldVals().size(), fieldVals());
                regIOobject::store(fieldVals);
            }
            break;

            case VTK_FLOAT:
            case VTK_DOUBLE:
            {
                auto fieldVals = autoPtr<scalarIOField>::New
                (
                    IOobject(arrayName, "", obj),
                    size
                );
                readBlock(inFile, fieldVals().size(), fieldVals());
                regIOobject::store(fieldVals);
            }
            break;

            case VTK_STRING:
            {
                DebugInfo
                    << "Reading strings:" << size << nl;

                auto fieldVals = autoPtr<stringIOList>::New
                (
                    IOobject(arrayName, "", obj),
                    size
                );
                // Consume current line.
                inFile.getLine(fieldVals()[0]);

                // Read without parsing
                for (string& s : fieldVals())
                {
                    inFile.getLine(s);
                }
                regIOobject::store(fieldVals);
            }
            break;

            default:
            {
                IOWarningInFunction(inFile)
                    << "Unhandled type " << dataType << nl
                    << "Skipping " << size
                    << " words." << nl;
                scalarField fieldVals;
                readBlock(inFile, size, fieldVals);
            }
            break;
        }
    }
    else
    {
        IOWarningInFunction(inFile)
            << "Unhandled type " << dataType << nl
            << "Skipping " << size
            << " words." << nl;
        scalarField fieldVals;
        readBlock(inFile, size, fieldVals);
    }
}


Foam::wordList Foam::vtkUnstructuredReader::readFieldArray
(
    ISstream& inFile,
    objectRegistry& obj,
    const label wantedSize
) const
{
    DynamicList<word> fields;

    word dataName(inFile);

    DebugInfo
        << "dataName:" << dataName << nl;

    label numArrays(readLabel(inFile));
    if (debug)
    {
        Pout<< "numArrays:" << numArrays << nl;
    }
    for (label i = 0; i < numArrays; i++)
    {
        word arrayName(inFile);
        label numComp(readLabel(inFile));
        label numTuples(readLabel(inFile));
        word dataType(inFile);

        DebugInfo
            << "Reading field " << arrayName
            << " of " << numTuples << " tuples of rank " << numComp << nl;

        if (wantedSize != -1 && numTuples != wantedSize)
        {
            FatalIOErrorInFunction(inFile)
                << "Expected " << wantedSize << " tuples but only have "
                << numTuples << nl
                << exit(FatalIOError);
        }

        readField
        (
            inFile,
            obj,
            arrayName,
            dataType,
            numTuples*numComp
        );
        fields.append(arrayName);
    }
    return fields.shrink();
}


Foam::objectRegistry& Foam::vtkUnstructuredReader::selectRegistry
(
    const parseMode readMode
)
{
    if (readMode == CELL_DATA)
    {
        return cellData_;
    }
    else if (readMode == POINT_DATA)
    {
        return pointData_;
    }

    return otherData_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkUnstructuredReader::vtkUnstructuredReader
(
    const objectRegistry& obr,
    ISstream& is
)
:
    version_(2.0),
    cellData_(IOobject("cellData", obr)),
    pointData_(IOobject("pointData", obr)),
    otherData_(IOobject("otherData", obr))
{
    read(is);
}


void Foam::vtkUnstructuredReader::read(ISstream& inFile)
{
    inFile.getLine(header_);
    DebugInfo<< "Header   : " << header_ << nl;

    // Extract version number from "# vtk DataFile Version 5.1"
    const auto split = stringOps::splitSpace(header_);
    if (split.size() >= 4 && split[split.size() - 2] == "Version")
    {
        float ver(0);
        if (readFloat(split[split.size() - 1], ver))
        {
            version_ = ver;
            DebugInfo<< "Version  : " << version_ << nl;
        }
    }

    inFile.getLine(title_);
    DebugInfo<< "Title    : " << title_ << nl;

    inFile.getLine(dataType_);
    DebugInfo<< "dataType : " << dataType_ << nl;

    if (dataType_ == "BINARY")
    {
        FatalIOErrorInFunction(inFile)
            << "Binary reading not supported" << nl
            << exit(FatalIOError);
    }

    parseMode readMode = NOMODE;
    label wantedSize = -1;


    // Intermediate storage for vertices of cells.
    labelList cellOffsets, cellVerts;

    token tok;

    while (inFile.read(tok).good() && tok.isWord())
    {
        const word tag = tok.wordToken();

        DebugInfo
            << "line:" << inFile.lineNumber()
            << " tag:" << tag << nl;

        if (tag == "DATASET")
        {
            word geomType(inFile);
            DebugInfo<< "geomType : " << geomType << nl;

            readMode = parseModeNames[geomType];
            wantedSize = -1;
        }
        else if (tag == "POINTS")
        {
            const label nPoints(readLabel(inFile));
            points_.resize(nPoints);

            DebugInfo
                << "Reading " << nPoints << " coordinates" << nl;

            word primitiveTag(inFile);
            if (primitiveTag != "float" && primitiveTag != "double")
            {
                FatalIOErrorInFunction(inFile)
                    << "Expected 'float' entry, found "
                    << primitiveTag << nl
                    << exit(FatalIOError);
            }
            for (point& p : points_)
            {
                inFile >> p.x() >> p.y() >> p.z();
            }
        }
        else if (tag == "CELLS")
        {
            label nCells(readLabel(inFile));
            const label nNumbers(readLabel(inFile));

            if (version_ < 5.0)
            {
                // VTK 4.2 and earlier (single-block)
                DebugInfo
                    << "Reading " << nCells
                    << " cells/faces (single block)" << nl;

                cellOffsets.clear();
                readBlock(inFile, nNumbers, cellVerts);
            }
            else
            {
                // VTK 5.0 and later (OFFSETS and CONNECTIVITY)

                const label nOffsets(nCells);
                --nCells;

                DebugInfo
                    << "Reading offsets/connectivity for "
                    << nCells << " cells/faces" << nl;

                readOffsetsConnectivity
                (
                    inFile,
                    "CELLS",
                    nOffsets, cellOffsets,
                    nNumbers, cellVerts
                );
            }
        }
        else if (tag == "CELL_TYPES")
        {
            const label nCellTypes(readLabel(inFile));

            labelList cellTypes;
            readBlock(inFile, nCellTypes, cellTypes);

            if (!cellTypes.empty() && cellVerts.empty())
            {
                FatalIOErrorInFunction(inFile)
                    << "Found " << cellTypes.size()
                    << " cellTypes but no cells." << nl
                    << exit(FatalIOError);
            }

            extractCells(inFile, cellTypes, cellOffsets, cellVerts);
            cellOffsets.clear();
            cellVerts.clear();
        }
        else if (tag == "LINES")
        {
            label nLines(readLabel(inFile));
            const label nNumbers(readLabel(inFile));

            labelList elemOffsets, elemVerts;

            if (version_ < 5.0)
            {
                // VTK 4.2 and earlier (single-block)
                DebugInfo
                    << "Reading " << nLines
                    << " lines (single block)" << nl;

                readBlock(inFile, nNumbers, elemVerts);
            }
            else
            {
                // VTK 5.0 and later (OFFSETS and CONNECTIVITY)

                const label nOffsets(nLines);
                --nLines;

                DebugInfo
                    << "Reading offsets/connectivity for "
                    << nLines << " lines" << nl;

                readOffsetsConnectivity
                (
                    inFile,
                    "LINES",
                    nOffsets, elemOffsets,
                    nNumbers, elemVerts
                );
            }

            // Append into lines
            label linei = lines_.size();
            lines_.resize(linei+nLines);
            lineMap_.resize(lines_.size());

            label dataIndex = 0;
            for (label i = 0; i < nLines; i++)
            {
                const label nVerts =
                (
                    elemOffsets.empty()
                  ? elemVerts[dataIndex++]
                  : (elemOffsets[i+1] - elemOffsets[i])
                );
                const labelSubList verts(elemVerts, nVerts, dataIndex);
                dataIndex += nVerts;

                lineMap_[linei] = linei;
                lines_[linei++] = verts;
            }
        }
        else if (tag == "POLYGONS")
        {
            label nFaces(readLabel(inFile));
            const label nNumbers(readLabel(inFile));

            labelList elemOffsets, elemVerts;

            if (version_ < 5.0)
            {
                // VTK 4.2 and earlier (single-block)
                DebugInfo
                    << "Reading " << nFaces
                    << " faces (single block)" << nl;

                readBlock(inFile, nNumbers, elemVerts);
            }
            else
            {
                // VTK 5.0 and later (OFFSETS and CONNECTIVITY)

                const label nOffsets(nFaces);
                --nFaces;

                DebugInfo
                    << "Reading offsets/connectivity for "
                    << nFaces << " faces" << nl;

                readOffsetsConnectivity
                (
                    inFile,
                    "POLYGONS",
                    nOffsets, elemOffsets,
                    nNumbers, elemVerts
                );
            }

            // Append into faces
            label facei = faces_.size();
            faces_.resize(facei+nFaces);
            faceMap_.resize(faces_.size());

            label dataIndex = 0;
            for (label i = 0; i < nFaces; ++i)
            {
                const label nVerts =
                (
                    elemOffsets.empty()
                  ? elemVerts[dataIndex++]
                  : (elemOffsets[i+1] - elemOffsets[i])
                );
                const labelSubList verts(elemVerts, nVerts, dataIndex);
                dataIndex += nVerts;

                faceMap_[facei] = facei;
                static_cast<labelList&>(faces_[facei++]) = verts;
            }
        }
        else if (tag == "POINT_DATA")
        {
            // 'POINT_DATA 24'
            readMode = POINT_DATA;
            wantedSize = points_.size();

            label nPoints(readLabel(inFile));
            if (nPoints != wantedSize)
            {
                FatalIOErrorInFunction(inFile)
                    << "Reading POINT_DATA : expected " << wantedSize
                    << " but read " << nPoints << nl
                    << exit(FatalIOError);
            }
        }
        else if (tag == "CELL_DATA")
        {
            readMode = CELL_DATA;
            wantedSize = cells_.size()+faces_.size()+lines_.size();

            label nCells(readLabel(inFile));
            if (nCells != wantedSize)
            {
                FatalIOErrorInFunction(inFile)
                    << "Reading CELL_DATA : expected "
                    << wantedSize
                    << " but read " << nCells << nl
                    << exit(FatalIOError);
            }
        }
        else if (tag == "FIELD")
        {
            // wantedSize already set according to type we expected to read.
            readFieldArray(inFile, selectRegistry(readMode), wantedSize);
        }
        else if (tag == "SCALARS")
        {
            string line;
            inFile.getLine(line);
            ISpanStream is(line);

            word dataName(is);
            word dataType(is);
            //label numComp(readLabel(inFile));

            DebugInfo
                << "Reading scalar " << dataName
                << " of type " << dataType
                << " from lookup table" << nl;

            word lookupTableTag(inFile);
            if (lookupTableTag != "LOOKUP_TABLE")
            {
                FatalIOErrorInFunction(inFile)
                    << "Expected tag LOOKUP_TABLE but read "
                    << lookupTableTag << nl
                    << exit(FatalIOError);
            }

            word lookupTableName(inFile);

            readField
            (
                inFile,
                selectRegistry(readMode),
                dataName,
                dataType,
                wantedSize//*numComp
            );
        }
        else if (tag == "VECTORS" || tag == "NORMALS")
        {
            // 'NORMALS Normals float'
            string line;
            inFile.getLine(line);
            ISpanStream is(line);

            word dataName(is);
            word dataType(is);
            DebugInfo
                << "Reading vector " << dataName
                << " of type " << dataType << nl;

            objectRegistry& reg = selectRegistry(readMode);

            readField
            (
                inFile,
                reg,
                dataName,
                dataType,
                3*wantedSize
            );

            if
            (
                vtkDataTypeNames[dataType] == VTK_FLOAT
             || vtkDataTypeNames[dataType] == VTK_DOUBLE
            )
            {
                objectRegistry::iterator iter = reg.find(dataName);
                scalarField s(*dynamic_cast<const scalarField*>(iter()));
                reg.erase(iter);
                auto fieldVals = autoPtr<vectorIOField>::New
                (
                    IOobject(dataName, "", reg),
                    s.size()/3
                );

                label elemI = 0;
                for (vector& val : fieldVals())
                {
                    val.x() = s[elemI++];
                    val.y() = s[elemI++];
                    val.z() = s[elemI++];
                }
                regIOobject::store(fieldVals);
            }
        }
        else if (tag == "TEXTURE_COORDINATES")
        {
            // 'TEXTURE_COORDINATES TCoords 2 float'
            string line;
            inFile.getLine(line);
            ISpanStream is(line);

            word dataName(is);          //"Tcoords"
            label dim(readLabel(is));
            word dataType(is);

            DebugInfo
                << "Reading texture coords " << dataName
                << " dimension " << dim
                << " of type " << dataType << nl;

            scalarField coords(dim*points_.size());
            readBlock(inFile, coords.size(), coords);
        }
        else if (tag == "TRIANGLE_STRIPS")
        {
            label nStrips(readLabel(inFile));
            const label nNumbers(readLabel(inFile));

            labelList elemOffsets, elemVerts;

            // Total number of faces in strips
            label nFaces = 0;

            if (version_ < 5.0)
            {
                // VTK 4.2 and earlier (single-block)
                DebugInfo
                    << "Reading " << nStrips
                    << " strips (single block)" << nl;

                readBlock(inFile, nNumbers, elemVerts);

                // Count number of faces (triangles)
                {
                    label dataIndex = 0;
                    for (label i = 0; i < nStrips; ++i)
                    {
                        const label nVerts = elemVerts[dataIndex++];
                        nFaces += nVerts-2;
                        dataIndex += nVerts;
                    }
                }
            }
            else
            {
                // VTK 5.0 and later (OFFSETS and CONNECTIVITY)

                const label nOffsets(nStrips);
                --nStrips;

                DebugInfo
                    << "Reading offsets/connectivity for "
                    << nStrips << " triangle strips." << nl;

                readOffsetsConnectivity
                (
                    inFile,
                    "TRIANGLE_STRIPS",
                    nOffsets, elemOffsets,
                    nNumbers, elemVerts
                );

                // Count number of faces (triangles)
                for (label i = 0; i < nStrips; ++i)
                {
                    const label nVerts = (elemOffsets[i+1] - elemOffsets[i]);
                    nFaces += nVerts-2;
                }
            }

            // Append into faces
            label facei = faces_.size();
            faces_.resize(facei+nFaces);
            faceMap_.resize(faces_.size());

            label dataIndex = 0;
            for (label i = 0; i < nStrips; ++i)
            {
                const label nVerts =
                (
                    elemOffsets.empty()
                  ? elemVerts[dataIndex++]
                  : (elemOffsets[i+1] - elemOffsets[i])
                );
                const label nTris = nVerts-2;

                // Advance before the first triangle
                if (nTris > 0)
                {
                    dataIndex += 2;
                }

                for (label triI = 0; triI < nTris; ++triI)
                {
                    faceMap_[facei] = facei;
                    face& f = faces_[facei++];

                    f.resize(3);

                    // NOTE: not clear if the orientation is correct
                    if ((triI % 2) == 0)
                    {
                        // Even (eg, 0-1-2, 2-3-4)
                        f[0] = elemVerts[dataIndex-2];
                        f[1] = elemVerts[dataIndex-1];
                    }
                    else
                    {
                        // Odd (eg, 2-1-3, 4-3-5)
                        f[0] = elemVerts[dataIndex-1];
                        f[1] = elemVerts[dataIndex-2];
                    }
                    f[2] = elemVerts[dataIndex++];
                }
            }
        }
        else if (tag == "METADATA")
        {
            word infoTag(inFile);
            if (infoTag != "INFORMATION")
            {
                FatalIOErrorInFunction(inFile)
                    << "Unsupported tag "
                    << infoTag << nl
                    << exit(FatalIOError);
            }
            label nInfo(readLabel(inFile));
            DebugInfo
                << "Ignoring " << nInfo << " metadata information." << nl;

            // Consume rest of line
            inFile.getLine(nullptr);
            for (label i = 0; i < 2*nInfo; i++)
            {
                inFile.getLine(nullptr);
            }
        }
        else
        {
            FatalIOErrorInFunction(inFile)
                << "Unsupported tag "
                << tag << nl
                << exit(FatalIOError);
        }
    }


    // There is some problem with orientation of prisms - the point
    // ordering seems to be different for some exports (e.g. of cgns)
    {
        const cellModel& prism = cellModel::ref(cellModel::PRISM);

        label nSwapped = 0;

        for (cellShape& shape : cells_)
        {
            if (shape.model() == prism)
            {
                const triPointRef bottom
                (
                    points_[shape[0]],
                    points_[shape[1]],
                    points_[shape[2]]
                );
                const triPointRef top
                (
                    points_[shape[3]],
                    points_[shape[4]],
                    points_[shape[5]]
                );

                const point bottomCc(bottom.centre());
                const vector bottomNormal(bottom.areaNormal());
                const point topCc(top.centre());

                if (((topCc - bottomCc) & bottomNormal) < 0)
                {
                    // Flip top and bottom
                    std::swap(shape[0], shape[3]);
                    std::swap(shape[1], shape[4]);
                    std::swap(shape[2], shape[5]);
                    ++nSwapped;
                }
            }
        }
        if (nSwapped)
        {
            WarningInFunction << "Swapped " << nSwapped
                << " prismatic cells" << nl;
        }
    }

    if (debug)
    {
        Info<< "Read points:" << points_.size()
            << " cells:" << cells_.size()
            << " faces:" << faces_.size()
            << " lines:" << lines_.size()
            << nl << nl;

        Info<< "Cell fields:" << nl;
        printFieldStats<vectorIOField>(cellData_);
        printFieldStats<scalarIOField>(cellData_);
        printFieldStats<labelIOField>(cellData_);
        printFieldStats<stringIOList>(cellData_);
        Info<< nl << nl;

        Info<< "Point fields:" << nl;
        printFieldStats<vectorIOField>(pointData_);
        printFieldStats<scalarIOField>(pointData_);
        printFieldStats<labelIOField>(pointData_);
        printFieldStats<stringIOList>(pointData_);
        Info<< nl << nl;

        Info<< "Other fields:" << nl;
        printFieldStats<vectorIOField>(otherData_);
        printFieldStats<scalarIOField>(otherData_);
        printFieldStats<labelIOField>(otherData_);
        printFieldStats<stringIOList>(otherData_);
    }
}


// ************************************************************************* //
