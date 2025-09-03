/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::vtk::polyWriter

Description
    Write faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkPolyWriter.C
    foamVtkPolyWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_polyWriter_H
#define Foam_vtk_polyWriter_H

#include "foamVtkFileWriter.H"
#include "pointField.H"
#include "edgeList.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class vtk::polyWriter Declaration
\*---------------------------------------------------------------------------*/

class polyWriter
:
    public vtk::fileWriter
{
protected:

    // Protected Data

        //- The number of field points for the current Piece
        label numberOfPoints_;

        //- The number of field cells (edges or faces) for the current Piece
        label numberOfCells_;

        //- Local number of points
        label nLocalPoints_;

        //- Local number of vertices (points)
        label nLocalVerts_;

        //- Local number of lines (edges)
        label nLocalLines_;

        //- Local number of polys (faces)
        label nLocalPolys_;


    // Protected Member Functions

        //- Write a uniform field of Cell (Poly or Line) or Point values
        template<class Type>
        void writeUniformValue
        (
            const label nCellValues,  // Could be Poly or Line!
            const word& fieldName,
            const Type& val
        );

private:

    // Private Member Functions

        //- Determine sizes (nLocalPoints_, nLocalLines_),
        //- and begin piece
        void beginPiece(const pointField& points, const edgeList& edges);

        //- Determine sizes (nLocalPoints_, nLocalPolys_),
        //- and begin piece
        void beginPiece(const pointField& points, const faceList& faces);

        //- Write points
        void writePoints(const pointField& points);

        //- Write lines, legacy format
        //  \param pointOffset processor-local point offset
        void writeLinesLegacy(const edgeList& edges, const label pointOffset);

        //- Write lines
        //  \param pointOffset processor-local point offset
        void writeLines(const edgeList& edges, const label pointOffset);

        //- Write faces, legacy format
        //  \param pointOffset processor-local point offset
        void writePolysLegacy(const faceList& faces, const label pointOffset);

        //- Write faces
        //  \param pointOffset processor-local point offset
        void writePolys(const faceList& faces, const label pointOffset);

        //- No copy construct
        polyWriter(const polyWriter&) = delete;

        //- No copy assignment
        void operator=(const polyWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        explicit polyWriter
        (
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        explicit polyWriter
        (
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        polyWriter
        (
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~polyWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::POLY_DATA);
        }

        //- Dummy write mesh topology method - Fatal if called
        //  Provided for special-purpose writers that do not hold
        //  their own geomety, but use writePolyGeometry() directly
        virtual bool writeGeometry();

        //- Low-level write edge/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writeLineGeometry
        (
            const pointField& points,
            const edgeList& edges
        );

        //- Low-level write face/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writePolyGeometry
        (
            const pointField& points,
            const faceList& faces
        );

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields = 0);

        //- Begin PointData for specified number of fields.
        //  Must be called prior to writing any point data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | CELL_DATA).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields = 0);


    // Write

        //- Write primitive field of CellData (Poly or Line) or PointData values
        template<class Type>
        void write(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of CellData
        template<class Type>
        void writeCellData(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of PointData
        template<class Type>
        void writePointData(const word& fieldName, const UList<Type>& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPolyWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
