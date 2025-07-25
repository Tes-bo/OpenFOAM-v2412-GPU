/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::vtk::coordSetWriter

Description
    Write as points/lines, optionally with fields,
    as a vtp file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Non-parallel only

SourceFiles
    foamVtkCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_coordSetWriter_H
#define Foam_vtk_coordSetWriter_H

#include "foamVtkPolyWriter.H"
#include "UPtrList.H"
#include "Field.H"
#include "instant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::coordSetWriter Declaration
\*---------------------------------------------------------------------------*/

class coordSetWriter
:
    public vtk::polyWriter
{
public:

    // Data Types

        //- Types of element output.
        enum elemOutputType
        {
            NO_ELEMENTS = 0,
            DEFAULT_ELEMENTS = 1,
            POINT_ELEMENTS = 2,
            LINE_ELEMENTS = 3
        };

private:

    // Private Member Data

        //- Reference to the points
        UPtrList<const pointField> points_;

        //- Time name/value
        instant instant_;

        //- Preferred output element type
        elemOutputType elemOutput_;


    // Private Member Functions


    // Private Member Functions

        //- Determine sizes (nLocalPoints_, nLocalLines_),
        //- and begin piece
        void beginPiece();

        //- Write points
        void writePoints();

        //- Write vertices, legacy format
        //  \param pointOffset processor-local point offset
        void writeVertsLegacy();

        //- Write vertices
        //  \param pointOffset processor-local point offset
        void writeVerts();

        //- Write lines, legacy format
        //  \param pointOffset processor-local point offset
        void writeLinesLegacy();

        //- Write lines
        //  \param pointOffset processor-local point offset
        void writeLines();


        //- No copy construct
        coordSetWriter(const coordSetWriter&) = delete;

        //- No copy assignment
        void operator=(const coordSetWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        coordSetWriter
        (
            const UPtrList<const pointField>& pts,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        coordSetWriter
        (
            const UPtrList<const pointField>& pts,
            const fileName& file,
            bool paralel = false  /* ignored: always false */
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        coordSetWriter
        (
            const UPtrList<const pointField>& pts,
            const vtk::outputOptions opts,
            const fileName& file,
            bool paralel = false  /* ignored: always false */
        );


    //- Destructor
    virtual ~coordSetWriter() = default;


    // Member Functions

        //- Define preferred element type
        void setElementType(const elemOutputType elemOutput);

        //- Open file for writing. Non-parallel only
        virtual bool open
        (
            const fileName& file,
            bool parallel = false  /* ignored: always false */
        );

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

        //- Reset point references to begin a new piece
        void piece(const UPtrList<const pointField>& points);


        //- Write processor ids for each line as CellData or for each point
        //- as PointData, depending on isPointData() state. No-op in serial.
        //  Not implemented.
        bool writeProcIDs();


    // Write

        //- Write primitive field of PointData
        template<class Type>
        void writePointData
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkCoordSetWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
