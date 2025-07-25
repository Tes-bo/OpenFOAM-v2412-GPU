/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::vtk::lagrangianWriter

Description
    Write lagrangian (cloud) positions and fields (as PointData) in
    VTP format. Legacy VTK format is intentionally not supported since
    the VTP format provides much better field selection in ParaView, and for
    consistency with the Foam::functionObjects::vtkCloud function object.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    If fields should be CellData instead of PointData (default), this
    must be decided at construction time.

SourceFiles
    lagrangianWriter.C
    lagrangianWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_lagrangianWriter_H
#define Foam_vtk_lagrangianWriter_H

#include "fvMesh.H"
#include "pointField.H"
#include "foamVtkFileWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;
template<class Type> class IOField;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class vtk::lagrangianWriter Declaration
\*---------------------------------------------------------------------------*/

class lagrangianWriter
:
    public vtk::fileWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- The cloud name
        const word cloudName_;

        //- The number of field points for the current Piece
        label numberOfPoints_;

        //- Write as CellData (verts) instead of as PointData.
        const bool useVerts_;


    // Private Member Functions

        //- The cloud directory for the current cloud name.
        fileName cloudDir() const;

        //- Transcribe the cloud into pointField
        pointField positions() const;

        //- Write vertex (cells)
        void writeVerts();


        //- No copy construct
        lagrangianWriter(const lagrangianWriter&) = delete;

        //- No copy assignment
        void operator=(const lagrangianWriter&) = delete;


protected:

    // Protected Member Functions

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is the number of fields, which is required for
        //      legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields=0);

        //- Begin PointData for specified number of fields.
        //  Must be called prior to writing any point data fields.
        //  \param nFields is the number of fields, which is required for
        //      legacy format.
        //  \note Expected calling states: (PIECE | CELL_DATA).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields=0);


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        //  \param useVerts Define VERTS and use CellData instead of PointData.
        lagrangianWriter
        (
            const fvMesh& mesh,
            const word& cloudName,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64,
            bool useVerts = false
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        lagrangianWriter
        (
            const fvMesh& mesh,
            const word& cloudName,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        lagrangianWriter
        (
            const fvMesh& mesh,
            const word& cloudName,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~lagrangianWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type. Always ".vtp"
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


        //- Begin parcel (PointData) output section
        //  Must be called prior to writing data fields.
        //  \note Expected calling states: (PIECE).
        //
        //  \return True if the state changed
        bool beginParcelData();

        //- Explicitly end parcel (PointData) output and switch to PIECE state
        //  Ignored (no-op) if not currently in the parcel state.
        bool endParcelData();


    // Write

        //- Write the IOField
        template<class Type>
        void write(const IOField<Type>& field);

        //- Write IOFields
        template<class Type>
        label writeFields(const wordList& fieldNames, bool verbose=true);

        //- Write IOFields
        template<class Type>
        label writeFields(const IOobjectList& objects, bool verbose=true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkLagrangianWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
