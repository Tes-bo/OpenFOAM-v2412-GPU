/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::faePatchField

Description
    faePatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation procedures, particularly with respect to specific fields.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faePatchField.C
    faePatchFieldBase.C
    faePatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faePatchField_H
#define Foam_faePatchField_H

#include "faPatch.H"
#include "DimensionedField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class objectRegistry;
class faPatchFieldMapper;
class edgeMesh;

template<class Type> class faePatchField;
template<class Type> class calculatedFaePatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faePatchField<Type>&);


/*---------------------------------------------------------------------------*\
                      Class faePatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for faePatchField
class faePatchFieldBase
{
    // Private Data

        //- Reference to patch
        const faPatch& patch_;

protected:

    // Protected Member Functions

        //- Read dictionary entries.
        //  Useful when initially constructed without a dictionary
        virtual void readDict(const dictionary& dict);

        //- Set updated state. This is a no-op for faePatchField
        void setUpdated(bool state) noexcept
        {}

        //- Set matrix manipulated state. This is a no-op for faePatchField
        void setManipulated(bool state) noexcept
        {}


public:

    //- Debug switch to disallow the use of generic faePatchField
    static int disallowGenericPatchField;

    //- Runtime type information
    TypeName("faePatchField");


    // Constructors

        //- Construct from patch
        explicit faePatchFieldBase(const faPatch& p);

        //- Construct from patch and patch type
        explicit faePatchFieldBase(const faPatch& p, const word& patchType);

        //- Construct from patch and dictionary
        faePatchFieldBase(const faPatch& p, const dictionary& dict);

        //- Copy construct with new patch
        faePatchFieldBase(const faePatchFieldBase& rhs, const faPatch& p);

        //- Copy construct
        faePatchFieldBase(const faePatchFieldBase& rhs);


    //- Destructor
    virtual ~faePatchFieldBase() = default;


    // Static Member Functions

        //- The type name for \c empty patch fields
        static const word& emptyType() noexcept
        {
            return Foam::fieldTypes::emptyType;
        }

        //- The type name for \c calculated patch fields
        static const word& calculatedType() noexcept
        {
            return Foam::fieldTypes::calculatedType;
        }


    // Member Functions

    // Attributes

        //- True if the patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        virtual bool fixesValue() const
        {
            return false;
        }

        //- True if the patch field is coupled
        virtual bool coupled() const
        {
            return false;
        }


    // Access

        //- The associated objectRegistry
        const objectRegistry& db() const;

        //- Return the patch
        const faPatch& patch() const noexcept
        {
            return patch_;
        }


    // Solution

        //- True if the boundary condition has already been updated.
        //- This is always true for faePatchField
        bool updated() const noexcept
        {
            return true;
        }

        //- True if the matrix has already been manipulated.
        //- Currently ignored (always false) for faePatchField
        bool manipulatedMatrix() const noexcept
        {
            return false;
        }


    // Check

        //- Check that patches are identical
        void checkPatch(const faePatchFieldBase& rhs) const;
};


/*---------------------------------------------------------------------------*\
                        Class faePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faePatchField
:
    public faePatchFieldBase,
    public Field<Type>
{
public:

    // Public Data Types

        //- The patch type for the patch field
        typedef faPatch Patch;

        //- The value_type for the patch field
        typedef Type value_type;

        //- The internal field type associated with the patch field
        typedef DimensionedField<Type, edgeMesh> Internal;

        //- Type for a \em calculated patch
        typedef calculatedFaePatchField<Type> Calculated;


private:

    // Private Data

        //- Reference to internal field
        const DimensionedField<Type, edgeMesh>& internalField_;

protected:

    // Protected Member Functions

        //- Read the "value" entry into \c *this.
        //  The reading can be optional (default), mandatory etc.
        //  \returns True on success
        bool readValueEntry
        (
            const dictionary& dict,
            IOobjectOption::readOption readOpt = IOobjectOption::LAZY_READ
        );

        //- Write \c *this field as a "value" entry
        void writeValueEntry(Ostream& os) const
        {
            Field<Type>::writeEntry("value", os);
        }


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            patch,
            (
                const faPatch& p,
                const DimensionedField<Type, edgeMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            patchMapper,
            (
                const faePatchField<Type>& ptf,
                const faPatch& p,
                const DimensionedField<Type, edgeMesh>& iF,
                const faPatchFieldMapper& m
            ),
            (dynamic_cast<const faePatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            dictionary,
            (
                const faPatch& p,
                const DimensionedField<Type, edgeMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch, internal field and value
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const Type& value
        );

        //- Construct from patch, internal field and patch field
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const Field<Type>& pfld
        );

        //- Construct from patch, internal field and patch field
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            Field<Type>&& pfld
        );

        //- Construct from patch, internal field and dictionary
        //  \note older versions have always treated "value" as optional
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary& dict,
            //! The "value" entry (default: mandatory)
            IOobjectOption::readOption requireValue = IOobjectOption::LAZY_READ
        );

        //- Construct by mapping the given faePatchField onto a new patch
        faePatchField
        (
            const faePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faePatchField(const faePatchField<Type>&);

        //- Construct as copy setting internal field reference
        faePatchField
        (
            const faePatchField<Type>&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Clone patch field with its own internal field reference
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>
            (
                new faePatchField<Type>(*this, this->internalField_)
            );
        }

        //- Clone with an internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, edgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type>>::New(*this, iF);
        }


    //- Destructor
    virtual ~faePatchField() = default;


    // Factory Methods

        //- Clone a patch field, optionally with internal field reference etc.
        template<class DerivedPatchField, class... Args>
        static tmp<faePatchField<Type>> Clone
        (
            const DerivedPatchField& pf,
            Args&&... args
        )
        {
            return tmp<faePatchField<Type>>
            (
                new DerivedPatchField(pf, std::forward<Args>(args)...)
            );
        }

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<faePatchField<Type>> New
        (
            const word& patchFieldType,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        //  Allows override of constraint type
        static tmp<faePatchField<Type>> New
        (
            const word& patchFieldType,
            const word& actualPatchType,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given faePatchField mapped onto a new patch
        static tmp<faePatchField<Type>> New
        (
            const faePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<faePatchField<Type>> New
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaePatchField created on
        //  freestore without setting patchField values
        static tmp<faePatchField<Type>> NewCalculatedType
        (
            const faPatch& p
        );

        //- Return a pointer to a new calculatedFaePatchField created on
        //  freestore without setting patchField values
        template<class AnyType>
        static tmp<faePatchField<Type>> NewCalculatedType
        (
            const faePatchField<AnyType>& pf
        );


    // Member Functions

    // Access

        //- Return const-reference to the dimensioned internal field
        const DimensionedField<Type, edgeMesh>& internalField() const noexcept
        {
            return internalField_;
        }

        //- Return const-reference to the internal field values
        const Field<Type>& primitiveField() const noexcept
        {
            return internalField_;
        }


        // Evaluation Functions

            //- Initialise the evaluation of the patch field after a local
            //  operation
            virtual void initEvaluateLocal
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            )
            {}

            //- Evaluate the patch field after a local operation (e.g. *=)
            virtual void evaluateLocal
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            )
            {}


        // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const faPatchFieldMapper&
        );

        //- Reverse map the given faePatchField onto this faePatchField
        virtual void rmap
        (
            const faePatchField<Type>&,
            const labelList&
        );


    // Evaluation Functions

        //- Initialise the evaluation of the patch field, generally a no-op
        virtual void initEvaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
        )
        {}

        //- Evaluate the patch field, generally a no-op
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
        )
        {}


    // Other

        //- Write
        virtual void write(Ostream& os) const;

        //- Check against given patch field
        void check(const faePatchField<Type>&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faePatchField<Type>&);
        virtual void operator+=(const faePatchField<Type>&);
        virtual void operator-=(const faePatchField<Type>&);
        virtual void operator*=(const faePatchField<scalar>&);
        virtual void operator/=(const faePatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const faePatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const faePatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;


    // Ostream Operator

        friend Ostream& operator<< <Type>(Ostream&, const faePatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faePatchField.C"
    #include "faePatchFieldNew.C"
    #include "calculatedFaePatchField.H"
#endif

// Runtime selection macros
#include "faePatchFieldMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
