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

Class
    Foam::vtk::internalWriter

Description
    Write an OpenFOAM volume (internal) geometry and internal fields
    as a vtu file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

See Also
    Foam::vtk::internalMeshWriter

SourceFiles
    foamVtkInternalWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_internalWriter_H
#define Foam_vtk_internalWriter_H

#include "foamVtkInternalMeshWriter.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class volPointInterpolation;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::internalWriter Declaration
\*---------------------------------------------------------------------------*/

class internalWriter
:
    public vtk::internalMeshWriter
{
    // Private Member Functions

        //- No copy construct
        internalWriter(const internalWriter&) = delete;

        //- No copy assignment
        void operator=(const internalWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64),
        internalWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        )
        :
            vtk::internalMeshWriter(mesh, cells, opts)
        {}


        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::internalMeshWriter(mesh, cells, file, parallel)
        {}

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::internalMeshWriter(mesh, cells, opts, file, parallel)
        {}


    //- Destructor
    virtual ~internalWriter() = default;


    // Member Functions

    // Write

        //- Write point field
        //  Interpolate to originating cell centre for decomposed cells.
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, pointMesh>& field
        );

        //- Write the internal field (CellData)
        template<class Type>
        void write
        (
            const DimensionedField<Type, volMesh>& field
        );

        //- Write the volume field (internal part)
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, volMesh>& field
        );

        //- Write internal field with point interpolation
        template<class Type>
        void write
        (
            const DimensionedField<Type, volMesh>& vfield,
            const volPointInterpolation& pInterp
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& vfield,
            const volPointInterpolation& pInterp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkInternalWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
