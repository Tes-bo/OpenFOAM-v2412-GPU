/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::vtk::GenericPatchGeoFieldsWriter

Description
    A Foam::vtk::GenericPatchWriter with support for geometric fields.

    Caution, currently only works properly for indirect patches.

SourceFiles
    foamVtkGenericPatchGeoFieldsWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_GenericPatchGeoFieldsWriter_H
#define Foam_vtk_GenericPatchGeoFieldsWriter_H

#include "areaFields.H"
#include "surfaceFields.H"
#include "foamVtkGenericPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
              Class vtk::GenericPatchGeoFieldsWriter Declaration
\*---------------------------------------------------------------------------*/

template<class PatchType>
class GenericPatchGeoFieldsWriter
:
    public vtk::GenericPatchWriter<PatchType>
{
    // Private Member Functions

        //- Get face field (internal face or boundary face)
        //- at specified face locations
        template<class Type>
        tmp<Field<Type>> getFaceField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
            const labelUList& faceAddr
        ) const;


    //- No copy construct
    GenericPatchGeoFieldsWriter
    (
        const GenericPatchGeoFieldsWriter<PatchType>&
    ) = delete;

    //- No copy assignment
    void operator=(const GenericPatchGeoFieldsWriter<PatchType>&) = delete;

public:

    // Constructors

        //- Inherit constructors from GenericPatchWriter
        using vtk::GenericPatchWriter<PatchType>::GenericPatchWriter;


    // Member Functions

        //- Write volume-mesh surface field (CellData) from specified mesh
        //- faces.
        //  For an indirect patch this will often simply be its addressing()
        template<class Type>
        void write
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
            const labelUList& faceAddr
        );

        //- Write area-mesh surface field (CellData)
        template<class Type>
        void write
        (
            const GeometricField<Type, faPatchField, areaMesh>& field
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkGenericPatchGeoFieldsWriter.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
