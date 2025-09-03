/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::vtk::surfaceWriter

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
    foamVtkSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_surfaceWriter_H
#define Foam_vtk_surfaceWriter_H

#include "foamVtkPolyWriter.H"
#include "instant.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceWriter
:
    public vtk::polyWriter
{
    // Private Member Data

        //- Reference to the points
        std::reference_wrapper<const pointField> points_;

        //- Reference to the faces
        std::reference_wrapper<const faceList> faces_;

        //- Time name/value
        instant instant_;


    // Private Member Functions

        //- No copy construct
        surfaceWriter(const surfaceWriter&) = delete;

        //- No copy assignment
        void operator=(const surfaceWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        surfaceWriter
        (
            const pointField& pts,
            const faceList& faces,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        surfaceWriter
        (
            const pointField& pts,
            const faceList& faces,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        surfaceWriter
        (
            const pointField& pts,
            const faceList& faces,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~surfaceWriter() = default;


    // Member Functions

        //- Define a time name/value for the output
        virtual void setTime(const instant& inst);

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "");

        //- Write patch topology
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();

        //- Write "TimeValue" FieldData (name as per Catalyst output)
        //  Must be called within the FIELD_DATA state.
        //  \note As a convenience this can also be called from
        //      (OPENED | DECLARED) states, in which case it invokes
        //      beginFieldData(1) internally.
        using vtk::fileWriter::writeTimeValue;

        //- Write the currently set time as "TimeValue" FieldData
        void writeTimeValue();

        //- Reset point/face references to begin a new piece
        void piece(const pointField& points, const faceList& faces);


        //- Write processor ids for each poly as CellData or for each point
        //- as PointData, depending on isPointData() state. No-op in serial.
        bool writeProcIDs();

        //- Write a uniform field of Cell (Poly) or Point values
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val)
        {
            polyWriter::writeUniformValue<Type>(nLocalPolys_, fieldName, val);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
