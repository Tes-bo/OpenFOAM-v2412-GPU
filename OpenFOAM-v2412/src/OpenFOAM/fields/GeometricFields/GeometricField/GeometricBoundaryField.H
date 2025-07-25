/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::GeometricBoundaryField

Description
    Generic GeometricBoundaryField class.

SourceFiles
    GeometricBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricBoundaryField_H
#define Foam_GeometricBoundaryField_H

#include "dimensionedTypes.H"
#include "DimensionedField.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class GeometricBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricBoundaryField
:
    public FieldField<PatchField, Type>
{
public:

    // Public Typedefs

        //- The boundary mesh type for the boundary fields
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- The internal field type associated with the boundary fields
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- This boundary field type
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> Boundary;

        //- The patch field type for the boundary fields
        typedef PatchField<Type> Patch;


private:

    // Private Data

        //- Reference to BoundaryMesh for which this field is defined
        const BoundaryMesh& bmesh_;


    // Private Member Functions

        //- Helper: check if patchfields have been evaluated. If not:
        //      exit = true  : FatalError
        //      exit = false : return bool
        template<class CheckPatchField>
        bool checkConsistency(const scalar tol, const bool exitIfBad) const;


public:

    //- Enable debug
    static int debug;

    //- User-defined tolerance (for consistency checks)
    static scalar tolerance;

    //- Enable local consistency
    static int localConsistency;


    // Constructors

        //- Construct from a BoundaryMesh, setting patches later
        explicit GeometricBoundaryField(const BoundaryMesh& bmesh);

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and a patch type
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and a wordList of patch types and optional the actual patch
        //- types (to override constraint patches)
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and a PtrList<PatchField<Type>> (to be cloned)
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Construct as copy, setting the reference to the internal field
        GeometricBoundaryField
        (
            const Internal& iField,
            const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf
        );

        //- Construct as copy, setting the reference to the internal field
        //- and resetting type of field for given patch IDs
        GeometricBoundaryField
        (
            const Internal& iField,
            const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf,
            const labelList& patchIDs,
            const word& patchFieldName
        );

        //- Copy construct
        //  Dangerous because Field may be set to a field which gets deleted
        //  Need new type of BoundaryField, one which is part of a geometric
        //  field for which snGrad etc. may be called and a free standing
        //  BoundaryField for which such operations are unavailable.
        GeometricBoundaryField(const GeometricBoundaryField& btf);

        //- Construct from dictionary
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const dictionary& dict
        );


    // Member Functions

        //- Read the boundary field
        void readField(const Internal& iField, const dictionary& dict);

        //- Update the boundary condition coefficients
        void updateCoeffs();

        //- Evaluate boundary conditions for each patch field.
        //- Uses specified or default comms.
        void evaluate
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions for patch fields matching the
        //- given predicate. Uses specified or default comms.
        template<class UnaryPredicate>
        void evaluate_if
        (
            const UnaryPredicate& pred,
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions after change in local values.
        //- Uses specified or default comms.
        void evaluateLocal
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions on coupled patches of given type.
        //- Uses specified or default comms.
        template<class CoupledPatchType>
        void evaluateCoupled
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions for selected patches.
        //- Uses non-blocking comms.
        void evaluateSelected(const labelUList& patchIDs);

        //- Return a list of the patch types
        wordList types() const;

        //- Return boundary field of values neighbouring the boundary
        tmp<GeometricBoundaryField> boundaryInternalField() const;

        //- Return a list of pointers for each patch field with only those
        //- pointing to interfaces being set
        LduInterfaceFieldPtrsList<Type> interfaces() const;

        //- Return a list of pointers for each patch field with only those
        //- pointing to interfaces being set
        lduInterfaceFieldPtrsList scalarInterfaces() const;

        //- Helper: check if field has been evaluated. See instantiations.
        bool check() const;


    // Write

        //- Write boundary field as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write dictionary entries of the individual boundary fields.
        void writeEntries(Ostream& os) const;


    // Member Operators

        //- Copy assignment from GeometricBoundaryField
        void operator=(const GeometricBoundaryField& bf);

        //- Copy assignment from FieldField\<PatchField, Type\>
        void operator=(const FieldField<PatchField, Type>& bf);

        //- Assignment to uniform value
        void operator=(const Type& val);

        //- Forced assignment from GeometricBoundaryField
        void operator==(const GeometricBoundaryField& bf);

        //- Forced assignment from FieldField\<PatchField, Type\>
        void operator==(const FieldField<PatchField, Type>& bf);

        //- Forced assignment to uniform value
        void operator==(const Type& val);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const GeometricBoundaryField&) = delete;
        bool operator!=(const FieldField<PatchField, Type>&) = delete;
        bool operator!=(const Type&) = delete;
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GeometricBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
