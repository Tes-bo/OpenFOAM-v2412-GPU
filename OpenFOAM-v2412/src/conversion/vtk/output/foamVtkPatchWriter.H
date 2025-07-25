/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::vtk::patchWriter

Description
    Write OpenFOAM patches and patch fields in VTP or legacy vtk format.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

See Also
    Foam::vtk::patchMeshWriter

SourceFiles
    foamVtkPatchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_patchWriter_H
#define Foam_vtk_patchWriter_H

#include "foamVtkPatchMeshWriter.H"
#include "volFields.H"
#include "pointFields.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class vtk::patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
:
    public vtk::patchMeshWriter
{
    // Private Member Data

        //- Use internal field value instead of patch value
        bool useNearCellValue_;


    // Private Member Functions

        //- No copy construct
        patchWriter(const patchWriter&) = delete;

        //- No copy assignment
        void operator=(const patchWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        //
        //  \param useNearCellValue to use cell instead of patch values
        patchWriter
        (
            const polyMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64,
            const bool useNearCellValue = false
        )
        :
            vtk::patchMeshWriter(mesh, patchIDs, opts),
            useNearCellValue_(useNearCellValue)
        {}

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        patchWriter
        (
            const polyMesh& mesh,
            const labelList& patchIDs,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::patchMeshWriter(mesh, patchIDs, file, parallel),
            useNearCellValue_(false)
        {}

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        patchWriter
        (
            const polyMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::patchMeshWriter(mesh, patchIDs, opts, file, parallel),
            useNearCellValue_(false)
        {}

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        //
        //  \param useNearCellValue to use cell instead of patch values
        patchWriter
        (
            const polyMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts,
            const bool useNearCellValue,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::patchMeshWriter(mesh, patchIDs, opts, file, parallel),
            useNearCellValue_(useNearCellValue)
        {}


    //- Destructor
    virtual ~patchWriter() = default;


    // Member Functions

    // Write

        //- Write point field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, pointMesh>& field
        );

        //- Write volume field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, volMesh>& field
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>& pInterp
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>* pInterp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPatchWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
