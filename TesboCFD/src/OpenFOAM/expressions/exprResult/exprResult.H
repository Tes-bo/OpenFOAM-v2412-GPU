/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::expressions::exprResult

Description
    A polymorphic field/result from evaluating an expression

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        resultType  | The type of result                    | no  | exprResult
        unsetValue  | Create without reading the dictionary | no  | false
        noReset     | Suppress reset on time                | no  | false
    \endtable

    When creating with values
    \table
        Property    | Description                           | Required | Default
        valueType   | Result value type (scalar, vector,..) | yes |
        isSingleValue | A single field value                | no  | false
        isPointValue  | Interpret values as point values    | no  | false
        value       | The field values                      | yes |
        fieldSize   | The size of the field (when not single-value) | no  |
    \endtable

SourceFiles
    exprResult.C
    exprResultI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResult_H
#define Foam_expressions_exprResult_H

#include "exprTraits.H"
#include "exprValue.H"
#include "dimensionedType.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprResult Declaration
\*---------------------------------------------------------------------------*/

class exprResult
{
    // Private Data

        //- The value type as string,
        //- normally corresponds to pTraits or typeName
        word valType_;

        //- Represents point data
        bool isPointData_;

        //- Whether or not the variable will be reset
        bool noReset_;

        //- Allow override of noReset_, but only accessible for subclasses
        bool needsReset_;

        //- Size (length) of field or object
        label size_;

        //- The single value representation
        expressions::exprValue value_;

        //- Allocated plain field (eg, scalarField)
        void *fieldPtr_;


    // Private Member Functions

        //- Type-checked deletion of the value pointer.
        //  \return True if the type check was satisfied
        template<class Type>
        inline bool deleteChecked();

        //- Invoke type-checked pointer deletion
        void destroy();

        //- Type-checked creation of field from dictionary (primitive) entry
        //  \return True if the type check was satisfied
        template<class Type>
        inline bool readChecked
        (
            const entry& e,
            //! Field length
            const label len,
            //! Read a single entry, not a field
            const bool singleValueOnly
        );

        //- Type-checked retrieval of uniform field from current results
        //  \return True if the type check was satisfied
        template<class Type>
        bool getUniformChecked
        (
            exprResult& result,
            const label size,
            const bool noWarn,
            const bool parRun
        ) const;

        //- Type-checked retrieval of \c bool uniform field from current result
        //  \return True if the type check was satisfied
        bool getUniformCheckedBool
        (
            exprResult& result,
            const label size,
            const bool noWarn,
            const bool parRun
        ) const;

        //- Type-checked determination of centre value (min/max)
        //  \return True if the type check was satisfied
        template<class Type>
        bool setAverageValueChecked(const bool parRun = UPstream::parRun());

        //- Type-checked determination of average bool value
        //  \return True if the type check was satisfied
        bool setAverageValueCheckedBool(const bool parRun = UPstream::parRun());

        //- Type-checked copy of field
        //  \return True if the type check was satisfied
        template<class Type>
        bool duplicateFieldChecked(const void* ptr);

        //- Type-checked writing of the single value (uniform) entry
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeSingleValueChecked(Ostream& os) const;

        //- Type-checked writing field as entry (if keyword is non-empty)
        //- or as plain field (if keyword is empty)
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeFieldChecked(const word& keyword, Ostream& os) const;

        //- Type-checked forwarding to Field::writeEntry
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeEntryChecked(const word& keyword, Ostream& os) const;

        //- Type-checked field addition with another expression field
        //  \return True if the type check was satisfied
        template<class Type>
        bool plusEqChecked(const exprResult& b);

        //- Type-checked field multiplication with a scalar
        //  \return True if the type check was satisfied
        template<class Type>
        bool multiplyEqChecked(const scalar& b);


        template<class Type>
        inline void setResultImpl(Field<Type>*, bool wantPointData=false);

        template<class Type>
        inline void setResultImpl(const Field<Type>&, bool wantPointData=false);

        template<class Type>
        inline void setResultImpl(Field<Type>&&, bool wantPointData=false);

        template<class Type>
        inline void setResultImpl(const Type& val, const label len);

        template<class Type>
        inline void setSingleValueImpl(const Type& val);


protected:

    // Protected Member Functions

        //- Simulate virtual templated methods
        inline virtual expressions::exprResult& target() { return *this; }

        //- Reset at new timestep according to the derived class type
        virtual void resetImpl();

        //- Reset at new timestep according to type
        //  \return true if it was actually reset
        bool reset(bool force=false);

        //- Adjusts the internal needsReset value
        void needsReset(bool val) { needsReset_ = val; }


public:

        //- An empty result
        static const exprResult null;

        //- Friendship with globals
        friend class exprResultGlobals;


        //- Runtime type information
        TypeName("exprResult");

        declareRunTimeSelectionTable
        (
            autoPtr,
            exprResult,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            exprResult,
            empty,
            (),
            ()
        );


    // Constructors

        //- Default construct
        exprResult();

        //- Copy construct
        exprResult(const exprResult& expr);

        //- Move construct
        exprResult(exprResult&& expr);

        //- Construct from a dictionary
        explicit exprResult
        (
            const dictionary& dict,
            //! Read a single entry, not a field
            const bool singleValueOnly = false,
            //! The "value" entry required or optional
            const bool valueReqd = false
        );

        //- Construct from Istream as dictionary content
        explicit exprResult(Istream& is);

        //- Construct by copying a field
        template<class Type>
        explicit exprResult(const Field<Type>& fld);

        //- Construct by moving a field
        template<class Type>
        explicit exprResult(Field<Type>&& fld);

        //- Construct for an IOobject
        template<class Type>
        explicit exprResult(autoPtr<Type>&& obj);

        //- Construct from a dimensioned value
        template<class Type>
        explicit exprResult(const dimensioned<Type>& dt);

        #undef exprResult_Construct
        #define exprResult_Construct(Type)                                \
            /*! \brief Construct from single value of Type */             \
            explicit exprResult(const Type& val) : exprResult()           \
            {                                                             \
                setSingleValue(val);                                      \
            }

        exprResult_Construct(bool);
        exprResult_Construct(scalar);
        exprResult_Construct(vector);
        exprResult_Construct(tensor);
        exprResult_Construct(symmTensor);
        exprResult_Construct(sphericalTensor);

        #undef exprResult_Construct


    // Selectors

        //- Return a reference to the selected value driver
        static autoPtr<exprResult> New(const dictionary& dict);

        //- Construct from Istream as dictionary content
        static autoPtr<exprResult> New(Istream& is);

        //- Clone
        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>::New(*this);
        }


    //- Destructor
    virtual ~exprResult();


    // Member Functions

    // Access

        //- Has a value?
        inline bool hasValue() const;

        //- Basic type for the field or single value
        inline const word& valueType() const noexcept;

        //- True if representing point data,
        //- or test for same value as wantPointData argument
        inline bool isPointData(const bool wantPointData=true) const;

        //- True if single, uniform value
        inline bool isUniform() const;

        //- True if valueType corresponds to the given Type
        template<class Type>
        inline bool isType() const;

        //- Return a single value when isUniform() is true,
        //- or Zero when it is non-uniform or if the type mismatches,
        //- which means that it can generally be considered as failsafe.
        template<class Type>
        inline Type getValue() const;

        //- Return a read pointer to the field data if the type matches,
        //- nullptr otherwise. Can generally be considered as failsafe.
        template<class Type>
        inline const Field<Type>* getField() const;

        //- True if valueType is a bool
        inline bool is_bool() const;

        //- The field or object size
        inline label size() const;

        //- The address of the field data content.
        //  Fatal for unknown types.
        //  Used, for example, for python integration
        const void* dataAddress() const;


    // Edit

        //- Clear (zero) the result
        void clear();

        //- Change reset behaviour
        void noReset() noexcept { noReset_ = true; }

        //- Change reset behaviour
        void allowReset() noexcept { noReset_ = false; }

        //- Test if field corresponds to a single-value and thus uniform.
        //  Uses field min/max to establish uniformity.
        //  Test afterwards with isUniform()
        void testIfSingleValue(const bool parRun = UPstream::parRun());


    // Set results

        //- Set result field, taking ownership of the pointer
        template<class Type>
        inline void setResult(Field<Type>*, bool wantPointData=false);

        //- Set result field, taking copy of the field contents
        template<class Type>
        inline void setResult(const Field<Type>&, bool wantPointData=false);

        //- Set result field, moving field contents
        template<class Type>
        inline void setResult(Field<Type>&&, bool wantPointData=false);

        //- Set uniform result field of given size
        template<class Type>
        inline void setResult(const Type& val, const label size);

        //- Set single-value uniform result
        template<class Type>
        inline void setSingleValue(const Type& val);


    // Access/Get results

        //- Return const reference to the field
        template<class Type>
        inline const Field<Type>& cref() const;

        //- Return non-const reference to the field
        template<class Type>
        inline Field<Type>& ref();

        //- Return non-const reference to the field, casting away constness
        template<class Type>
        inline Field<Type>& constCast() const;

        //- Return tmp field of the contents,
        //- optionally keeping a copy in cache
        template<class Type>
        inline tmp<Field<Type>> getResult(bool cacheCopy=false);

        //- Construct a uniform field from the current results
        //  Uses the field average. Optionally warning if the min/max
        //  deviation is larger than SMALL.
        exprResult getUniform
        (
            const label size,
            const bool noWarn,
            const bool parRun = UPstream::parRun()
        ) const;

        //- Get a reduced result
        template<template<class> class BinaryOp, class Type>
        inline Type getReduced
        (
            const BinaryOp<Type>& bop,
            const Type& initial = pTraits<Type>::zero
        );


    // Write

        //- Forwarding to Field::writeEntry
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write entry as dictionary contents
        void writeDict(Ostream& os, const bool subDict=true) const;

        //- Write the field, optionally as an entry
        void writeField(Ostream& os, const word& keyword = "") const;

        //- Write the single value, or the first value from field
        void writeValue(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        virtual void operator=(const exprResult& rhs);

        //- Move assignment
        virtual void operator=(exprResult&& rhs);


        //- Scalar multiplication
        exprResult& operator*=(const scalar& b);

        //- Addition of results
        exprResult& operator+=(const exprResult& b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Operators

expressions::exprResult operator*
(
    const scalar& a,
    const expressions::exprResult& b
);
expressions::exprResult operator*
(
    const expressions::exprResult& a,
    const scalar& b
);
expressions::exprResult operator+
(
    const expressions::exprResult& a,
    const expressions::exprResult& b
);


// IO Operator
Istream& operator>>(Istream& os, expressions::exprResult& data);
Ostream& operator<<(Ostream& os, const expressions::exprResult& data);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exprResultI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
