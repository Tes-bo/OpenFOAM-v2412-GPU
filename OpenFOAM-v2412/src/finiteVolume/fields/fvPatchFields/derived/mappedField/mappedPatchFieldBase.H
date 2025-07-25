/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::mappedPatchFieldBase

Description
    Functionality for sampling fields using mappedPatchBase. Every call to
    mappedField() returns a sampled field, optionally scaled to maintain an
    area-weighted average.

    Example usage:
    \verbatim
    {
        field               T;          // default is same as fvPatchField
        setAverage          false;
        average             1.0;        // only if setAverage=true
        interpolationScheme cellPoint;  // default is cell
    }
    \endverbatim

    \table
        Property     | Description                      | Required | Default
        field        | name of field to be mapped       | no  | this field name
        setAverage   | use average value                | no  | false |
        average      | average value to apply if \c setAverage = yes | partly |
        interpolationScheme | interpolation scheme      | partly | cell |
    \endtable

SourceFiles
    mappedPatchFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchFieldBase_H
#define mappedPatchFieldBase_H

#include "fixedValueFvPatchFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class> class interpolation;
class mapDistribute;
class mappedPatchBase;

/*---------------------------------------------------------------------------*\
                    Class mappedPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedPatchFieldBase
{
    // Private Member Functions

        //- Selective retrieval of "average" entry from the dictionary
        static Type getAverage(const dictionary& dict, bool mandatory);


protected:

    // Protected Data

        //- Mapping engine
        const mappedPatchBase& mapper_;

        //- Underlying patch field
        const fvPatchField<Type>& patchField_;

        //- Name of field to sample
        word fieldName_;

        //- If true adjust the mapped field to maintain average value average_
        const bool setAverage_;

        //- Average value the mapped field is adjusted to maintain if
        //- setAverage_ is set true
        const Type average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;


    // Protected Member Functions

        //- Construct field from registered elements
        template<class T>
        void initRetrieveField
        (
            const objectRegistry& obr,
            const word& region,
            const word& patch,
            const labelListList& map,
            const word& fieldName,
            const Field<T>& fld
        ) const;

        //- Helper : storeField and retrieveField and interpolate. Leaves fld
        //  unchanged (and returns false) if new values cannot be retrieved.
        //  Returns true otherwise.
        template<class T>
        bool storeAndRetrieveField
        (
            const word& fieldName,
            const label myComm,
            const labelListList& subMap,
            const label constructSize,
            const labelListList& constructMap,
            const labelListList& address,
            const scalarListList& weights,
            Field<T>& fld
        ) const;


        //- Initialise field to retrieve (e.g. when value read from dictionary)
        template<class T>
        void initRetrieveField
        (
            const word& fieldName,
            const Field<T>& fld
        ) const;


        //- Store elements of field onto (sub) registry
        template<class T>
        void storeField
        (
            const objectRegistry& obr,
            const word& region,
            const word& patch,
            const label myComm,
            const labelListList& procToMap,
            const word& fieldName,
            const Field<T>& fld
        ) const;

        //- Construct field from registered elements.
        //  \return true if successful
        template<class T>
        bool retrieveField
        (
            const bool allowUnset,
            const objectRegistry& obr,
            const word& region,
            const word& patch,
            const label myComm,
            const labelListList& procToMap,
            const word& fieldName,
            Field<T>& fld
        ) const;


public:

    // Constructors

        //- Construct from components
        mappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct from dictionary
        mappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const dictionary& dict
        );

        //- Construct from dictionary and initial field value
        //  (for when running in database mode)
        mappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const dictionary& dict,
            const Field<Type>& fld
        );

        //- Construct empty
        mappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField
        );

        //- Construct copy
        explicit mappedPatchFieldBase
        (
            const mappedPatchFieldBase<Type>& mapper
        );

        //- Construct copy, resetting patch and field
        mappedPatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const mappedPatchFieldBase<Type>& base
        );


    //- Destructor
    virtual ~mappedPatchFieldBase() = default;


    // Member Functions

        //- Field to sample. Either on my or nbr mesh
        template<class T>
        const GeometricField<T, fvPatchField, volMesh>& sampleField
        (
            const word& fieldName
        ) const;

        //- Field to sample. Either on my or nbr mesh
        const GeometricField<Type, fvPatchField, volMesh>& sampleField() const;

        //- Map sampleField onto *this patch
        virtual tmp<Field<Type>> mappedField() const;

        //- Map internal of sampleField onto *this patch
        virtual tmp<Field<Type>> mappedInternalField() const;

        //- Map optional weightField onto *this patch. Default is deltaCoeffs
        virtual tmp<scalarField> mappedWeightField() const;

        //- Map optional weightField (given by name) onto *this patch
        virtual void mappedWeightField
        (
            const word& weightFieldName,
            tmp<scalarField>& thisWeights,
            tmp<scalarField>& nbrWeights
        ) const;

        //- Write
        virtual void write(Ostream& os) const;


    // Helpers

        //- Check that patch is of correct type
        static const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        );

        //- Wrapper for mapDistribute::distribute that knows about
        //- dabase mapping
        template<class T>
        void distribute(const word& fieldName, Field<T>& newValues) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchFieldBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
