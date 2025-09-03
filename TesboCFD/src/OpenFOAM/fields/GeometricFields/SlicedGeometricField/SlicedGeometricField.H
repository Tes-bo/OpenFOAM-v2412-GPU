/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::SlicedGeometricField

Description
    Specialization of GeometricField which holds slices of given complete
    fields in a form that they act as a GeometricField.

    The destructor is wrapped to avoid deallocation of the storage of the
    complete fields when this is destroyed.

    SlicedGeometricField can only be instantiated with a valid form of
    SlicedPatchField to handle the slicing and storage deallocation of the
    boundary field.

SourceFiles
    SlicedGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SlicedGeometricField_H
#define Foam_SlicedGeometricField_H

#include "GeometricField.H"
#include "SlicedDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SlicedGeometricField Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField
:
    public GeometricField<Type, PatchField, GeoMesh>
{
public:

    // Public Typedefs

        typedef typename GeoMesh::Mesh Mesh;
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- The internal field slice for this SlicedGeometricField
        typedef SlicedDimensionedField<Type, GeoMesh> Internal;


private:

    // Private Member Functions

        //- Return true if field size is too small to be the complete field
        //- (ie, it must be boundary field addressing)
        static bool
        isBoundaryAddressing(const Mesh& mesh, const label fieldSize);

        //- Slice the given field and a create a PtrList of SlicedPatchField
        //  from which the boundary field is built
        tmp<FieldField<PatchField, Type>>
        makeBoundary
        (
            const Mesh& mesh,
            const Field<Type>& completeOrBoundaryField,
            const bool preserveCouples,
            const bool preserveProcessorOnly = false,
            const bool isBoundaryOnly = false
        ) const;

        //- Slice the given field and a create a PtrList of SlicedPatchField
        //  from which the boundary field is built
        tmp<FieldField<PatchField, Type>>
        makeBoundary
        (
            const Mesh& mesh,
            const FieldField<PatchField, Type>& bField,
            const bool preserveCouples
        ) const;

        // Note - copy construct allowed

        //- No copy assignment
        void operator=(const SlicedGeometricField&) = delete;

        //- No copy assignment from GeometricField,
        //  == assignment is allowed.
        void operator=
        (
            const GeometricField<Type, PatchField, GeoMesh>&
        ) = delete;

        //- No copy assignment from tmp<GeometricField>,
        //  == assignment is allowed.
        void operator=
        (
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        ) = delete;


public:

    // Constructors

        //- Construct from components and field to slice
        SlicedGeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet& dims,
            const Field<Type>& completeField,
            const bool preserveCouples = true
        );

        //- Construct from components and separate fields to slice for the
        //- internal field and boundary field
        SlicedGeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet& dims,
            const Field<Type>& completeIField,
            const Field<Type>& completeBField,
            const bool preserveCouples = true,
            const bool preserveProcessorOnly = false
        );

        //- Construct from GeometricField. Reuses full internal and
        //- patch fields except on couples (preserveCouples=true).
        SlicedGeometricField
        (
            const IOobject&,
            const GeometricField<Type, PatchField, GeoMesh>&,
            const bool preserveCouples = true
        );

        //- Copy construct
        SlicedGeometricField
        (
            const SlicedGeometricField
            <
                Type,
                PatchField,
                SlicedPatchField,
                GeoMesh
            >&
        );

        //- Clone
        tmp<SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>>
        clone() const;


    //- Destructor
    ~SlicedGeometricField();


    // Member Functions

        //- Correct boundary field
        void correctBoundaryConditions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SlicedGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
