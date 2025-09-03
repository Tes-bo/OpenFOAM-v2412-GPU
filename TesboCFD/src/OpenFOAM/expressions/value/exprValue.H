/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::expressions::exprValue

Description
    A polymorphic typed union of simple primitive and VectorSpace types.
    It uses a 'fatter' representation that includes standard VectorSpace
    types, which avoids heap allocations at the expense of more storage.
    This is mostly not an issue since lists and fields would box/unbox
    an entire field, not individual values.

SourceFiles
    exprValue.C
    exprValueI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprValue_H
#define Foam_expressions_exprValue_H

#include "exprTraits.H"  // For valueTypeCode and label, scalar, vector etc.
#include "contiguous.H"
#include "InfoProxy.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Execute macro for known exprValue types, with more frequent ones first
#undef  FOR_ALL_EXPR_VALUE_TYPES
#define FOR_ALL_EXPR_VALUE_TYPES(Macro, ...)                                  \
    Macro(scalar, __VA_ARGS__)                                                \
    Macro(vector, __VA_ARGS__)                                                \
    Macro(tensor, __VA_ARGS__)                                                \
    Macro(symmTensor, __VA_ARGS__)                                            \
    Macro(sphericalTensor, __VA_ARGS__)                                       \
    Macro(label, __VA_ARGS__)                                                 \
    Macro(bool, __VA_ARGS__)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ITstream;

namespace expressions
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                       Class exprValueUnion Declaration
\*---------------------------------------------------------------------------*/

//- The data content (as a union).
//  The stem of the member names correspond identically to the data types
//  to enable macro definitions (naming similar to Foam::token)
union exprValueUnion
{
    // Member Data

        #undef  declareUnionMember
        #define declareUnionMember(Type, UnusedParam)                         \
                                                                              \
        /*! Define Type data member */                                        \
        Type Type##Value;

        FOR_ALL_EXPR_VALUE_TYPES(declareUnionMember);
        #undef declareUnionMember


    // Member Functions

        //- Runtime 'assert' for unimplemented generic methods
        //  This seems to be the best way to catch programming errors
        //  since static_assert does not help here.
        //  The 'noexcept' is slightly misleading (needed for the caller)
        static void notSpecialized(const std::string& msg) noexcept;

        //- Return read pointer to typed union member,
        //- which is nullptr for unspecialized versions
        template<class Type>
        inline const Type* get() const noexcept { return nullptr; }

        //- No generic set(), only specialized methods
        template<class Type>
        inline void set(const Type& val) noexcept
        {
            notSpecialized("set<" + std::string(typeid(Type).name()) + '>');
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Type-specific get/set (definitions)
#undef  defineMultiTypeValueUnionMethods
#define defineMultiTypeValueUnionMethods(Type, UnusedParam)                   \
                                                                              \
/*! Return read pointer to Type union member */                               \
template<>                                                                    \
inline const Type* exprValueUnion::get<Type>() const noexcept                 \
{                                                                             \
    return &(Type##Value);                                                    \
}                                                                             \
                                                                              \
/*! Set Type value */                                                         \
template<>                                                                    \
inline void exprValueUnion::set<Type>(const Type& val) noexcept               \
{                                                                             \
    Type##Value = val;                                                        \
}

// Type-specific methods
FOR_ALL_EXPR_VALUE_TYPES(defineMultiTypeValueUnionMethods);
#undef defineMultiTypeValueUnionMethods

} // End namespace Detail

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                       Class exprValue Declaration
\*---------------------------------------------------------------------------*/

class exprValue
{
    // Private Data

        //- The data content (as a union)
        Detail::exprValueUnion data_;

        //- The data type (as a type-code)
        expressions::valueTypeCode typeCode_ = valueTypeCode::NONE;


    // Private Member Functions

        //- Fill the data with zero, preserving the data type
        inline void fill_zero();

        //- Copy assignment
        void deepCopy(const exprValue& rhs);

public:

    // Constructors

        //- Default construct (zero-initialized) as 'none'
        inline exprValue();

        //- Copy construct (deep copy)
        inline exprValue(const exprValue& rhs);

        //- Construct from Type. Fatal for unsupported types
        template<class Type>
        inline explicit exprValue(const Type& val);


    // Static Methods (valueTypeCode)

        //- True if valueTypeCode is not none/invalid
        inline static bool good
        (
            const expressions::valueTypeCode
        ) noexcept;

        //- True if valueTypeCode is bool/label
        inline static bool is_integral
        (
            const expressions::valueTypeCode
        ) noexcept;

        //- The number of components associated with the valueTypeCode
        inline static direction nComponents
        (
            const expressions::valueTypeCode
        ) noexcept;

        //- The vector-space rank associated with the valueTypeCode
        inline static direction rank
        (
            const expressions::valueTypeCode
        ) noexcept;

        //- True if the specified type is supported
        template<class Type>
        inline static bool supportedType();


    // Static Methods

        //- Detect the type from the available tokens.
        //  Possible ambiguity between label and scalar.
        //  \return INVALID for unsupported type
        static expressions::valueTypeCode peekType(const ITstream& is);

        //- Read entire string as a exprValue,
        //- skipping leading/trailing whitespace.
        //  \return True if successful.
        static bool read(const std::string& str, exprValue& val);


    // Member Functions

        //- The value type code
        expressions::valueTypeCode typeCode() const noexcept
        {
            return typeCode_;
        }

        //- True if the value type is not none/invalid
        inline bool good() const noexcept;

        //- True if the value type is an integral (bool or label)
        inline bool is_integral() const noexcept;

        //- The number of components associated with the value type
        inline direction nComponents() const noexcept;

        //- The vector-space rank associated with the value type
        inline direction rank() const noexcept;

        //- The name for the value type. Similar to pTraits typeName
        inline word valueTypeName() const;

        //- The name for the compound token (for a List of values).
        //- Eg, \c List\<scalar\>, \c List\<vector\> ...
        inline word listCompoundName() const;

        //- Reset to 'none'
        void clear();

        //- Read read tokens (if possible)
        //  \return True on success
        bool read(Istream& is);

        //- Guess type and read tokens (if possible)
        //  \return True on success
        bool readTokens(ITstream& is);

        //- Compare (type,value)
        int compare(const exprValue& rhs) const;


    // Typed Access

        //- Return non-null pointer to the data element (if types match).
        //- Can also be tested as a bool.
        template<class Type>
        inline const Type* isA() const noexcept;

        //- Assign from type.
        //- Returns false and sets to 'none' for unsupported types
        template<class Type>
        inline bool set(const Type& val);

        //- Return the value. Return pTraits zero if the types mismatch
        template<class Type>
        inline const Type& get() const;


    // Type-specific methods, operators

    #undef  defineUnionMethods
    #define defineUnionMethods(Type, UnusedParam)                             \
                                                                              \
        /*! True if value is a Type */                                        \
        bool is_##Type() const noexcept                                       \
        {                                                                     \
            return (typeCode_ == exprTypeTraits<Type>::value);                \
        }                                                                     \
                                                                              \
        /*! Assign from Type */                                               \
        void operator=(const Type& val) { this->set<Type>(val); }

        FOR_ALL_EXPR_VALUE_TYPES(defineUnionMethods)
        #undef defineUnionMethods


        //- Copy assignment
        void operator=(const exprValue& rhs) { deepCopy(rhs); }

        //- Assign from zero. Changes value but not type
        void operator=(const Foam::zero) { fill_zero(); }


    // Low-level access

        //- Pointer to the data content as byte data
        const char* cdata_bytes() const noexcept
        {
            return reinterpret_cast<const char*>(this);
        }

        //- Pointer to the data content as byte data
        char* data_bytes() noexcept
        {
            return reinterpret_cast<char*>(this);
        }

        //- Size of the (contiguous) data content as byte data
        //- is compile-time constant
        static constexpr unsigned size_bytes() noexcept
        {
            return sizeof(exprValue);
        }


    // Member Functions

        //- Compare (type,value) for equality
        bool operator==(const exprValue& rhs) const;

        //- Compare (type,value) for inequality
        bool operator!=(const exprValue& rhs) const { return !(*this == rhs); }

        //- Compare (type,value)
        bool operator<(const exprValue& rhs) const;


    // Output

        //- Return info proxy for printing information to a stream
        InfoProxy<exprValue> info() const { return *this; }

        //- Write the (type-specific) content.
        //  For none : emits pair of brackets.
        //  For invalid : emits "bad".
        //
        //  Use OTstream for the stream to recover as tokens.
        //
        //  \param prune suppress the output for none/invalid
        void write(Ostream& os, bool prune = false) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions

//- The data content are contiguous
template<> struct is_contiguous<expressions::exprValue> : std::true_type {};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read/parse value from input stream (uses ASCII format).
Istream& operator>>(Istream& is, expressions::exprValue& val);

//- Write value to output stream (uses ASCII format).
//- Writes 'none' or 'bad' for unknown/unsupported types.
Ostream& operator<<(Ostream& os, const expressions::exprValue& val);

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<expressions::exprValue>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exprValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
