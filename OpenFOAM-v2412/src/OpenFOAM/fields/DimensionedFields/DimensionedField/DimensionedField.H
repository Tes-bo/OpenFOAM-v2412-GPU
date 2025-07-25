/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::DimensionedField

Description
    Field with dimensions and associated with geometry type GeoMesh which is
    used to size the field and a reference to it is maintained.

SourceFiles
    DimensionedFieldI.H
    DimensionedField.C
    DimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_DimensionedField_H
#define Foam_DimensionedField_H

#include "regIOobject.H"
#include "Field.H"
#include "dimensionedType.H"
#include "orientedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type, class GeoMesh> class DimensionedField;

template<class Type, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const DimensionedField<Type, GeoMesh>& df
);

template<class Type, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
);


/*---------------------------------------------------------------------------*\
                      Class DimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedField
:
    public regIOobject,
    public Field<Type>
{
public:

    // Public Typedefs

        //- The mesh type for the DimensionedField
        typedef typename GeoMesh::Mesh Mesh;

        //- The DimensionedField is usually synonymous with an internal field
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- Type of the field from which this DimensionedField is derived
        typedef Field<Type> FieldType;

        //- Component type of the field elements
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //- Reference to mesh
        const Mesh& mesh_;

        //- Dimension set for this field
        dimensionSet dimensions_;

        //- Oriented flag
        orientedType oriented_;


    // Private Member Functions

        //- Assert that non-zero field size == mesh size
        void checkFieldSize() const;

        //- Read from file if it is present. Checks the readOption flag
        void readIfPresent(const word& fieldDictEntry = "value");

        //- Read the field - create the field dictionary on-the-fly
        void readField(const word& fieldDictEntry = "value");

        //- Implementation for 'New' with specified registerObject preference.
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        template<class... Args>
        static tmp<DimensionedField<Type, GeoMesh>> New_impl
        (
            IOobjectOption::registerOption regOpt,
            const word& name,
            const Mesh& mesh,
            Args&&... args
        );

public:

    //- Runtime type information
    TypeName("DimensionedField");


    // Static Member Functions

        //- Return a null DimensionedField (reference to a nullObject).
        static const DimensionedField<Type, GeoMesh>& null() noexcept
        {
            return NullObjectRef<DimensionedField<Type, GeoMesh>>();
        }


    // Constructors

        //- Construct from components, copy initial field content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& field
        );

        //- Construct from components, transferring initial field content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& field
        );

        //- Construct from components, transferring initial field content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            List<Type>&& field
        );

        //- Construct from components, copy or transfer tmp content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const tmp<Field<Type>>& tfield
        );

        //- Construct from components, setting initial size and dimensions,
        //- but not initialising any field values.
        //  Used for temporary fields which are initialised after construction
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const bool checkIOFlags = true
        );

        //- Construct from components, setting dimensions and initial
        //- uniform field value.
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims,
            const bool checkIOFlags = true
        );

        //- Construct from components, setting dimensions and initial
        //- uniform field value.
        //  The name of the dimensioned\<Type\> has no influence.
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const bool checkIOFlags = true
        );

        //- Construct from Istream.
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const word& fieldDictEntry = "value"
        );

        //- Construct from dictionary
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dictionary& fieldDict,
            const word& fieldDictEntry = "value"
        );

        //- Copy construct
        DimensionedField(const DimensionedField<Type, GeoMesh>& df);

        //- Move construct
        DimensionedField(DimensionedField<Type, GeoMesh>&& df);

        //- Copy construct or reuse (move) as specified.
        DimensionedField(DimensionedField<Type, GeoMesh>& df, bool reuse);

        //- Construct from tmp\<DimensionedField\> deleting argument
        DimensionedField
        (
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Copy construct, resetting IO parameters
        DimensionedField
        (
            const IOobject& io,
            const DimensionedField<Type, GeoMesh>& df
        );

        //- Move construct, resetting IO parameters
        DimensionedField
        (
            const IOobject& io,
            DimensionedField<Type, GeoMesh>&& df
        );

        //- Copy or move construct, resetting IO parameters.
        DimensionedField
        (
            const IOobject& io,
            DimensionedField<Type, GeoMesh>& df,
            bool reuse
        );

        //- Construct from tmp\<DimensionedField\> deleting argument,
        //- resetting IO parameters.
        DimensionedField
        (
            const IOobject& io,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Copy construct with a new name
        DimensionedField
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh>& df
        );

        //- Move construct with a new name
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh>&& df
        );

        //- Copy or move construct, resetting name.
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh>& df,
            bool reuse
        );

        //- Construct with a new name from tmp\<DimensionedField\>
        DimensionedField
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Clone
        tmp<DimensionedField<Type, GeoMesh>> clone() const;


    // Factory Methods

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions, copy of internal field.
        //- [Takes current timeName from the mesh registry].
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            IOobjectOption::registerOption regOpt,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions, copy of internal field.
        //- [Takes current timeName from the mesh registry].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions, move internal field contents.
        //- [Takes current timeName from the mesh registry].
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            IOobjectOption::registerOption regOpt,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& iField
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions, move internal field contents.
        //- [Takes current timeName from the mesh registry].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& iField
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions.
        //- [Takes current timeName from the mesh registry].
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            IOobjectOption::registerOption regOpt,
            const Mesh& mesh,
            const dimensionSet& dims
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensions.
        //- [Takes current timeName from the mesh registry].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, value, dimensions.
        //- [Takes current timeName from the mesh registry].
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            IOobjectOption::registerOption regOpt,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, value, dimensions.
        //- [Takes current timeName from the mesh registry].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensioned-type.
        //- [Takes current timeName from the mesh registry].
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            IOobjectOption::registerOption regOpt,
            const Mesh& mesh,
            const dimensioned<Type>& dt
        );

        //- Return tmp field (NO_READ, NO_WRITE)
        //- from name, mesh, dimensioned-type.
        //- [Takes current timeName from the mesh registry].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>& dt
        );

        //- Return renamed tmp field (NO_READ, NO_WRITE)
        //- retaining its instance/local.
        //  For LEGACY_REGISTER, registration is determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            IOobjectOption::registerOption regOpt,
            const tmp<DimensionedField<Type, GeoMesh>>& tfld
        );

        //- Return renamed tmp field (NO_READ, NO_WRITE)
        //- retaining its instance/local.
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh>>& tfld
        );

        //- Construct tmp field (NO_READ, NO_WRITE)
        //- based on mesh/registry information from an existing field.
        //- [Takes instance from the field].
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        template<class AnyType>
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const DimensionedField<AnyType, GeoMesh>& fld,
            const word& name,
            const dimensionSet& dims
        );

        //- Construct tmp field (NO_READ, NO_WRITE)
        //- based on mesh/registry information from an existing field
        //- and initialise with value.
        //  Registration/persistence determined by
        //  objectRegistry::is_cacheTemporaryObject().
        template<class AnyType>
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const DimensionedField<AnyType, GeoMesh>& fld,
            const word& name,
            const dimensioned<Type>& dt
        );


    //- Destructor
    virtual ~DimensionedField();


    // Member Functions

        //- Resize to the mesh size and read the field from the dictionary.
        //- Internal implementation checks field vs mesh size.
        void readField
        (
            const dictionary& fieldDict,
            const word& fieldDictEntry = "value"
        );

        //- Return mesh
        inline const Mesh& mesh() const noexcept;

        //- Return dimensions
        inline const dimensionSet& dimensions() const noexcept;

        //- Return non-const access to dimensions
        inline dimensionSet& dimensions() noexcept;

        //- Return oriented type
        inline orientedType oriented() const noexcept;

        //- Return non-const access to the oriented type
        inline orientedType& oriented() noexcept;

        //- True if field is ORIENTED
        inline bool is_oriented() const noexcept;

        //- Set the oriented flag: on/off
        inline void setOriented(bool on = true) noexcept;

        //- Return const-reference to the field values
        inline const Field<Type>& field() const noexcept;

        //- Return reference to the field values
        inline Field<Type>& field() noexcept;

        //- Return a component field of the field
        tmp<DimensionedField<cmptType, GeoMesh>> component
        (
            const direction d
        ) const;

        //- Replace a component field of the field
        void replace
        (
            const direction d,
            const DimensionedField<cmptType, GeoMesh>& df
        );

        //- Replace a component field of the field
        void replace
        (
            const direction d,
            const tmp<DimensionedField<cmptType, GeoMesh>>& tdf
        );

        // Inherits clamp, clamp_min, clamp_max (without dimensions) from Field

        //- Return the field transpose (only defined for second rank tensors)
        tmp<DimensionedField<Type, GeoMesh>> T() const;

        //- Calculate and return arithmetic average
        dimensioned<Type> average() const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const DimensionedField<scalar, GeoMesh>& weightField
        ) const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const tmp<DimensionedField<scalar, GeoMesh>>& tweightField
        ) const;


    // Write

        //- Write dimensions, oriented flag (if valid) and the
        //- field data as a dictionary entry with the specified name.
        bool writeData(Ostream& os, const word& fieldDictEntry) const;

        //- The writeData function (required by regIOobject),
        //- calls writeData with dictionary entry name = "value"
        bool writeData(Ostream& os) const { return writeData(os, "value"); }


    // Member Operators

        void operator=(const DimensionedField<Type, GeoMesh>& df);
        void operator=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        //- Assign dimensions and value.
        void operator=(const dimensioned<Type>& dt);

        void operator+=(const DimensionedField<Type, GeoMesh>& df);
        void operator+=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        void operator-=(const DimensionedField<Type, GeoMesh>& df);
        void operator-=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        void operator*=(const DimensionedField<scalar, GeoMesh>& df);
        void operator*=(const tmp<DimensionedField<scalar, GeoMesh>>& tdf);

        void operator/=(const DimensionedField<scalar, GeoMesh>& df);
        void operator/=(const tmp<DimensionedField<scalar, GeoMesh>>& tdf);

        void operator+=(const dimensioned<Type>& dt);
        void operator-=(const dimensioned<Type>& dt);

        void operator*=(const dimensioned<scalar>& dt);
        void operator/=(const dimensioned<scalar>& dt);


    // Ostream Operators

        //- Calls DimensionedField::writeData()
        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream& os,
            const DimensionedField<Type, GeoMesh>& df
        );

        //- Calls DimensionedField::writeData()
        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream& os,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldI.H"
#include "DimensionedFieldFunctions.H"

#ifdef NoRepository
    #include "DimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
