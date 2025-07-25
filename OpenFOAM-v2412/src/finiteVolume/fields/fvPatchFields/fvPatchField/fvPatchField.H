/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fvPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation procedures, particularly with respect to specific
    fields.

SourceFiles
    fvPatchField.C
    fvPatchFieldBase.C
    fvPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvPatchField_H
#define Foam_fvPatchField_H

#include "fvPatch.H"
#include "DimensionedField.H"
#include "fieldTypes.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class objectRegistry;
class fvPatchFieldMapper;
class volMesh;

template<class Type> class fvPatchField;
template<class Type> class calculatedFvPatchField;
template<class Type> class fvMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                      Class fvPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for fvPatchField
class fvPatchFieldBase
{
    // Private Data

        //- Reference to patch
        const fvPatch& patch_;

        //- Update index used so that updateCoeffs is called only once during
        //- the construction of the matrix
        bool updated_;

        //- Update index used so that manipulateMatrix is called only once
        //- during the construction of the matrix
        bool manipulatedMatrix_;

        //- Use implicit formulation
        bool useImplicit_;

        //- Optional patch type
        //  Used to allow specified boundary conditions to be applied
        //  to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;


protected:

    // Protected Member Functions

        //- Read dictionary entries.
        //  Useful when initially constructed without a dictionary
        virtual void readDict(const dictionary& dict);


public:

        //- Set updated state
        void setUpdated(bool state) noexcept
        {
            updated_ = state;
        }

        //- Set matrix manipulated state
        void setManipulated(bool state) noexcept
        {
            manipulatedMatrix_ = state;
        }


public:

    //- Debug switch to disallow the use of generic fvPatchField
    static int disallowGenericPatchField;

    //- Runtime type information
    TypeName("fvPatchField");


    // Constructors

        //- Construct from patch
        explicit fvPatchFieldBase(const fvPatch& p);

        //- Construct from patch and patch type
        explicit fvPatchFieldBase(const fvPatch& p, const word& patchType);

        //- Construct from patch and dictionary
        fvPatchFieldBase(const fvPatch& p, const dictionary& dict);

        //- Copy construct with new patch
        fvPatchFieldBase(const fvPatchFieldBase& rhs, const fvPatch& p);

        //- Copy construct
        fvPatchFieldBase(const fvPatchFieldBase& rhs);


    //- Destructor
    virtual ~fvPatchFieldBase() = default;


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

        //- The type name for \c extrapolatedCalculated patch fields
        //- combines \c zero-gradient and \c calculated
        static const word& extrapolatedCalculatedType() noexcept
        {
            return Foam::fieldTypes::extrapolatedCalculatedType;
        }

        //- The type name for \c zeroGradient patch fields
        static const word& zeroGradientType() noexcept
        {
            return Foam::fieldTypes::zeroGradientType;
        }


    // Member Functions

    // Attributes

        //- True if the value of the patch field is altered by assignment
        virtual bool assignable() const
        {
            return true;
        }

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
        const fvPatch& patch() const noexcept
        {
            return patch_;
        }

        //- The optional patch type
        const word& patchType() const noexcept
        {
            return patchType_;
        }

        //- The optional patch type
        word& patchType() noexcept
        {
            return patchType_;
        }


    // Solution

        //- True if the boundary condition has already been updated
        bool updated() const noexcept
        {
            return updated_;
        }

        //- True if the matrix has already been manipulated
        bool manipulatedMatrix() const noexcept
        {
            return manipulatedMatrix_;
        }

        //- Use implicit formulation for coupled patches only
        bool useImplicit() const noexcept
        {
            return useImplicit_;
        }

        //- Set useImplicit on/off
        //  \return old value
        bool useImplicit(bool on) noexcept
        {
            bool old(useImplicit_);
            useImplicit_ = on;
            return old;
        }


    // Check

        //- Check that patches are identical
        void checkPatch(const fvPatchFieldBase& rhs) const;
};


/*---------------------------------------------------------------------------*\
                        Class fvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvPatchField
:
    public fvPatchFieldBase,
    public Field<Type>
{
public:

    // Public Data Types

        //- The patch type for the patch field
        typedef fvPatch Patch;

        //- The value_type for the patch field
        typedef Type value_type;

        //- The internal field type associated with the patch field
        typedef DimensionedField<Type, volMesh> Internal;

        //- Type for a \em calculated patch
        typedef calculatedFvPatchField<Type> Calculated;


private:

    // Private Data

        //- Reference to internal field
        const DimensionedField<Type, volMesh>& internalField_;


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

        //- Assign the patch field from the internal field
        void extrapolateInternal();


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and patch type
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const word& patchType
        );

        //- Construct from patch, internal field and value
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Type& value
        );

        //- Construct from patch, internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& pfld
        );

        //- Construct from patch, internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            Field<Type>&& pfld
        );

        //- Construct from patch, internal field and dictionary
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary& dict,
            //! The "value" entry (default: mandatory)
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Construct, forwarding to readOption variant
        fvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const dictionary& dict,
            const bool needValue
        )
        :
            fvPatchField
            (
                p, iF, dict,
                (needValue? IOobjectOption::MUST_READ : IOobjectOption::NO_READ)
            )
        {}

        //- Construct by mapping the given fvPatchField onto a new patch
        fvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvPatchField(const fvPatchField<Type>&);

        //- Construct as copy setting internal field reference
        fvPatchField
        (
            const fvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Clone patch field with its own internal field reference
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fvPatchField<Type>(*this, this->internalField_)
            );
        }

        //- Clone patch field with given internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fvPatchField<Type>(*this, iF)
            );
        }


    // Factory Methods

        //- Clone a patch field, optionally with internal field reference etc.
        template<class DerivedPatchField, class... Args>
        static tmp<fvPatchField<Type>> Clone
        (
            const DerivedPatchField& pf,
            Args&&... args
        )
        {
            return tmp<fvPatchField<Type>>
            (
                new DerivedPatchField(pf, std::forward<Args>(args)...)
            );
        }

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvPatchField<Type>> New
        (
            const word& patchFieldType,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<fvPatchField<Type>> New
        (
            const word& patchFieldType,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvPatchField mapped onto a new patch
        static tmp<fvPatchField<Type>> New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatch& p
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class AnyType>
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatchField<AnyType>& pf
        );


    //- Destructor
    virtual ~fvPatchField() = default;


    // Member Functions

    // Access

        //- Return const-reference to the dimensioned internal field
        const DimensionedField<Type, volMesh>& internalField() const noexcept
        {
            return internalField_;
        }

        //- Return const-reference to the internal field values
        const Field<Type>& primitiveField() const noexcept
        {
            return internalField_;
        }


        // Mapping Functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation Functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const;

            //- Return patch-normal gradient for coupled-patches
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            //  with a weight field (0..1). This weight field is usually
            //  provided as the amount of geometric overlap for 'duplicate'
            //  patches. Sets Updated to true
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Return internal field next to patch
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Extract internal field next to patch
            // \param [out] pfld The extracted patch field.
            virtual void patchInternalField(Field<Type>& pfld) const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField() const
            {
                NotImplemented;
                return *this;
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            )
            {}

            //- Evaluate the patch field, sets updated() to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            );

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

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
                return *this;
            }


            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& matrix,
                const scalarField& weights
            );

            //- Manipulate fvMatrix
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& matrix,
                const label iMatrix,
                const direction cmp
            );


    // Other

        //- Write
        virtual void write(Ostream&) const;

        //- Check against given patch field
        void check(const fvPatchField<Type>&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

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

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const fvPatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;


    // Ostream Operator

        friend Ostream& operator<< <Type>(Ostream&, const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvPatchField.C"
    #include "fvPatchFieldNew.C"
    #include "calculatedFvPatchField.H"
    #include "zeroGradientFvPatchField.H"
#endif

// Runtime selection macros
#include "fvPatchFieldMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
