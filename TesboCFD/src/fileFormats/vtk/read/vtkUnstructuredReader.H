/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Class
    Foam::vtkUnstructuredReader

Description
    Reader for vtk UNSTRUCTURED_GRID legacy files.
    Supports single CELLS, POINTS etc. entry only.

    - all integer types (int, unsigned_int, long etc.) become Foam::label
    - all real types (float, double) become Foam::scalar
    - POINTS becomes OpenFOAM points
    - CELLS gets split into OpenFOAM
        - cells
        - faces
        - lines
    - CELL_DATA or POINT_DATA gets stored on the corresponding objectRegistry
      in original vtk numbering order so use e.g. faceMap() to go from entry
      in faces() back to vtk numbering.

SourceFiles
    vtkUnstructuredReader.C

\*---------------------------------------------------------------------------*/

#ifndef vtkUnstructuredReader_H
#define vtkUnstructuredReader_H

#include "foamVtkCore.H"
#include "objectRegistry.H"
#include "cellShapeList.H"
#include "HashSet.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class vtkUnstructuredReader Declaration
\*---------------------------------------------------------------------------*/

class vtkUnstructuredReader
{
public:

    // Public Data Types

        //- Enumeration defining the vtk data types
        enum vtkDataType
        {
            VTK_INT,
            VTK_INT64,
            VTK_UINT,
            VTK_LONG,
            VTK_ULONG,
            VTK_FLOAT,
            VTK_DOUBLE,
            VTK_STRING,
            VTK_ID
        };

        static const Enum<vtkDataType> vtkDataTypeNames;


        //- Enumeration defining the vtk dataset types
        enum vtkDataSetType
        {
            VTK_FIELD,
            VTK_SCALARS,
            VTK_VECTORS
        };

        static const Enum<vtkDataSetType> vtkDataSetTypeNames;


        //- Enumeration defining the parse mode - type of data being read
        enum parseMode
        {
            NOMODE,
            UNSTRUCTURED_GRID,
            POLYDATA,
            CELL_DATA,
            POINT_DATA
        };

        static const Enum<parseMode> parseModeNames;


private:

    //- The VTK version
    float version_;

    //- Header
    string header_;

    //- Title
    string title_;

    //- DataType
    string dataType_;


    // Geometry

        //- Points
        pointField points_;

        //- 3D cells.
        cellShapeList cells_;

        //- Map from cells back to original ID
        labelList cellMap_;

        //- 2D cells (=faces)
        faceList faces_;

        //- Map from faces back to original ID
        labelList faceMap_;

        //- 1D cells (=edges)
        labelListList lines_;

        labelList lineMap_;


    // Data

        //- Cell based fields
        objectRegistry cellData_;

        //- Point based fields
        objectRegistry pointData_;

        //- Other fields
        objectRegistry otherData_;



    // Private Member Functions

        //- Read OFFSETS, CONNECTIVITY arrays
        static void readOffsetsConnectivity
        (
            ISstream& is,
            const char* entryName,
            const label nOffsets,
            labelList& offsets,
            const label nConnectivity,
            labelList& connectivity
        );

        static void warnUnhandledType
        (
            const Istream& is,  // For error message
            const label type,
            labelHashSet& warningGiven
        );

        //- Split cellTypes into cells, faces and lines
        void extractCells
        (
            const Istream& is,  // For error message
            const labelUList& cellTypes,
            const labelUList& cellOffsets,
            const labelUList& cellVertData
        );

        //- Read single field and stores it on the objectRegistry.
        void readField
        (
            ISstream& inFile,
            objectRegistry& obj,
            const word& arrayName,
            const word& dataType,
            const label size
        ) const;

        //- Reads fields, stores them on the objectRegistry. Returns a list of
        //  read fields
        wordList readFieldArray
        (
            ISstream& inFile,
            objectRegistry& obj,
            const label wantedSize
        ) const;

        objectRegistry& selectRegistry(const parseMode readMode);

        void read(ISstream& inFile);

        //- No copy assignment
        void operator=(const vtkUnstructuredReader&) = delete;


public:

    //- Runtime type information
    ClassName("vtkUnstructuredReader");


    // Constructors

        //- Construct from input stream, read all
        vtkUnstructuredReader(const objectRegistry& obr, ISstream& is);


    // Member Functions

        //- Header
        const string& header() const noexcept
        {
            return header_;
        }

        //- Title
        const string& title() const noexcept
        {
            return title_;
        }

        //- DataType
        const string& dataType() const noexcept
        {
            return dataType_;
        }


        //- Points
        const pointField& points() const noexcept
        {
            return points_;
        }

        pointField& points() noexcept
        {
            return points_;
        }

        //- 3D cells
        const cellShapeList& cells() const noexcept
        {
            return cells_;
        }

        cellShapeList& cells() noexcept
        {
            return cells_;
        }

        const labelList& cellMap() const noexcept
        {
            return cellMap_;
        }

        //- 2D cells (=faces)
        const faceList& faces() const noexcept
        {
            return faces_;
        }

        faceList& faces() noexcept
        {
            return faces_;
        }

        const labelList& faceMap() const noexcept
        {
            return faceMap_;
        }

        //- 1D cells (=open lines)
        const labelListList& lines() const noexcept
        {
            return lines_;
        }

        labelListList& lines() noexcept
        {
            return lines_;
        }

        const labelList& lineMap() const noexcept
        {
            return lineMap_;
        }

        //- Cell based fields
        const objectRegistry& cellData() const noexcept
        {
            return cellData_;
        }

        objectRegistry& cellData() noexcept
        {
            return cellData_;
        }

        //- Point based fields
        const objectRegistry& pointData() const noexcept
        {
            return pointData_;
        }

        objectRegistry& pointData() noexcept
        {
            return pointData_;
        }

        //- Other fields
        const objectRegistry& otherData() const noexcept
        {
            return otherData_;
        }

        objectRegistry& otherData() noexcept
        {
            return otherData_;
        }


        //- Debug: print contents of objectRegistry
        template<class Type>
        static void printFieldStats(const objectRegistry&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkUnstructuredReaderTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
