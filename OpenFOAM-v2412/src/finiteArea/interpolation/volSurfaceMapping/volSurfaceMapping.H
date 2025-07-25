/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    volSurfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_volSurfaceMapping_H
#define Foam_volSurfaceMapping_H

#include "areaFields.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "PtrList.H"
#include "PtrMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class fvPatchField;
template<class Type> class fvsPatchField;

/*---------------------------------------------------------------------------*\
                      Class volSurfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class volSurfaceMapping
{
    // Private Data

        //- Reference to mesh
        const faMesh& mesh_;


public:

    // Constructors

        //- Construct from mesh
        explicit volSurfaceMapping(const faMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- No copy construct
        volSurfaceMapping(const volSurfaceMapping&) = delete;

        //- No copy assignment
        void operator=(const volSurfaceMapping&) = delete;


    //- Destructor
    ~volSurfaceMapping() = default;


    // Mapping to area fields

        //- Map volume boundary fields as area field
        template<class Type>
        void mapToSurface
        (
            const GeometricBoundaryField<Type, fvPatchField, volMesh>&,
            Field<Type>& result
        ) const;

        //- Map volume boundary fields as area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const GeometricBoundaryField<Type, fvPatchField, volMesh>&
        ) const;

        //- Map volume (boundary) fields to area field
        template<class Type>
        void mapToSurface
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            Field<Type>& result
        ) const;

        //- Map volume (boundary) fields to area field
        template<class Type>
        void mapToSurface
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            Field<Type>& result
        ) const;

        //- Map volume (boundary) fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Map volume (boundary) fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;


        //- Map surface fields to area field
        template<class Type>
        void mapToSurface
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            Field<Type>& result
        ) const;

        //- Map surface fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;

        //- Map surface fields to area field
        template<class Type>
        void mapToSurface
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&,
            Field<Type>& result
        ) const;

        //- Map surface fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
        ) const;


        //- Map pre-calculated boundary fields to area field
        template<class Type>
        void mapToSurface
        (
            const UPtrList<Field<Type>>& patchFields,
            Field<Type>& result
        ) const;

        //- Map pre-calculated boundary fields to area field
        template<class Type>
        void mapToSurface
        (
            const PtrMap<Field<Type>>& patchFields,
            Field<Type>& result
        ) const;

        //- Map pre-calculated boundary fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const UPtrList<Field<Type>>& patchFields
        ) const;

        //- Map pre-calculated boundary fields to area field
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const PtrMap<Field<Type>>& patchFields
        ) const;


        //- Map patch internal field to area field
        template<class Type>
        void mapInternalToSurface
        (
            const GeometricBoundaryField<Type, fvPatchField, volMesh>&,
            Field<Type>& result
        ) const;

        //- Map patch internal field to area field
        template<class Type>
        tmp<Field<Type>> mapInternalToSurface
        (
            const GeometricBoundaryField<Type, fvPatchField, volMesh>&
        ) const;

        //- Map patch internal field to area field
        template<class Type>
        void mapInternalToSurface
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            Field<Type>& result
        ) const;

        //- Map patch internal field to area field
        template<class Type>
        void mapInternalToSurface
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            Field<Type>& result
        ) const;

        //- Map patch internal field to area field
        template<class Type>
        tmp<Field<Type>> mapInternalToSurface
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Map patch internal field to area field
        template<class Type>
        tmp<Field<Type>> mapInternalToSurface
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;


    // Mapping from area field to volume (boundary) field

        //- Map area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const DimensionedField<Type, areaMesh>&,
            GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;

        //- Map area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const tmp<DimensionedField<Type, areaMesh>>&,
            GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;

        //- Map area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const DimensionedField<Type, areaMesh>&,
            GeometricField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;

        //- Map area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const tmp<DimensionedField<Type, areaMesh>>&,
            GeometricField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;

        //- Map tmp area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
            GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;

        //- Map tmp area field to volume boundary field,
        //- optionally restricted to a single destination patch
        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
            GeometricField<Type, fvPatchField, volMesh>& dest,
            const label destPatchi = -1
        ) const;


    // Mapping from area field to volume patch field

        //- Map area field to a volume boundary patch
        template<class Type>
        void mapToVolumePatch
        (
            const DimensionedField<Type, areaMesh>& af,
            Field<Type>& dest,
            const label destPatchi
        ) const;

        //- Map tmp area field to a volume boundary patch
        template<class Type>
        void mapToVolumePatch
        (
            const tmp<DimensionedField<Type, areaMesh>>& af,
            Field<Type>& dest,
            const label destPatchi
        ) const;

        //- Map tmp area field to a volume boundary patch
        template<class Type>
        void mapToVolumePatch
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
            Field<Type>& dest,
            const label destPatchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volSurfaceMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
