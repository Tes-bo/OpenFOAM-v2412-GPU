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
    Foam::vtk::GenericPatchWriter

Description
    Write concrete PrimitivePatch faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.
    The patch type is limited to representations of polygon faces
    and 3D points and must support the following methods:

    - localPoints()
    - localFaces()
    .

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_GenericPatchWriter_H
#define Foam_vtk_GenericPatchWriter_H

#include "foamVtkPolyWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::GenericPatchWriter Declaration
\*---------------------------------------------------------------------------*/

template<class PatchType>
class GenericPatchWriter
:
    public vtk::polyWriter
{
    // Private Member Data

        //- Reference to faces/points as a patch
        const PatchType& pp_;


    // Private Member Functions

        //- No copy construct
        GenericPatchWriter(const GenericPatchWriter<PatchType>&) = delete;

        //- No copy assignment
        void operator=(const GenericPatchWriter<PatchType>&) = delete;


public:

    // Constructors

        //- Construct from patch (default output INLINE_BASE64)
        explicit GenericPatchWriter
        (
            const PatchType& pp,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        )
        :
            vtk::polyWriter(opts),
            pp_(pp)
        {}

        //- Construct from components (default output INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        GenericPatchWriter
        (
            const PatchType& pp,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::polyWriter(file, parallel),
            pp_(pp)
        {}

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        GenericPatchWriter
        (
            const PatchType& pp,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        )
        :
            vtk::polyWriter(opts, file, parallel),
            pp_(pp)
        {}


    //- Destructor
    virtual ~GenericPatchWriter() = default;


    // Member Functions

        //- Reference to the originating face/points patch
        const PatchType& patch() const noexcept
        {
            return pp_;
        }

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "surface")
        {
            return vtk::polyWriter::beginFile(title);
        }

        //- Write patch topology
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry()
        {
            return writePolyGeometry(pp_.localPoints(), pp_.localFaces());
        }


    // Write Fields

        //- Write processor ids for each poly as CellData
        bool writeProcIDs()
        {
            if (this->isPointData())
            {
                return vtk::fileWriter::writeProcIDs(nLocalPoints_);
            }
            return vtk::polyWriter::writeProcIDs(nLocalPolys_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
