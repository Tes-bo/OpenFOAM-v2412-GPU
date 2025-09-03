/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::vtk::surfaceFieldWriter

Description
    Write surfaces fields (as PointData) in VTP format.
    Legacy VTK format is intentionally not supported.

    The file output is structured as HEAD, FIELD_DATA, PIECE followed by any
    CELL_DATA or POINT_DATA. These states are tracked internally to help
    detect logic errors.

    The FieldData element is to be placed prior to writing any geometry
    Piece. This moves the information to the front of the output file
    for visibility and simplifies the logic.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    surfaceFieldWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_surfaceFieldWriter_H
#define Foam_vtk_surfaceFieldWriter_H

#include "foamVtkFileWriter.H"
#include "fvMesh.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class vtk::surfaceFieldWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceFieldWriter
:
    public vtk::fileWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- The number of field points for the current Piece
        label numberOfPoints_;


public:

    // Generated Methods

        //- No copy construct
        surfaceFieldWriter(const surfaceFieldWriter&) = delete;

        //- No copy assignment
        void operator=(const surfaceFieldWriter&) = delete;


    // Constructors

        //- Construct from mesh (default format INLINE_BASE64)
        surfaceFieldWriter
        (
            const fvMesh& mesh,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from mesh (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        surfaceFieldWriter
        (
            const fvMesh& mesh,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from mesh and open the file for writing.
        //  The file name is with/without an extension.
        surfaceFieldWriter
        (
            const fvMesh& mesh,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~surfaceFieldWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        inline static word ext(vtk::outputOptions)
        {
            // No legacy
            return vtk::fileExtension[vtk::fileTag::POLY_DATA];
        }


        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "");

        //- Write cloud positions
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is the number of fields, which is required for
        //      legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields=0);

        //- Begin PointData output section
        //  Must be called prior to writing data fields.
        //  \note Expected calling states: (PIECE).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields=0);


    // Write

        //- Write field
        void write(const surfaceVectorField& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
