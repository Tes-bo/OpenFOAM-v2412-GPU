/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::exprTypeTraits

Description
    Simple type identifiers for polymorphic expression values.
    The definitions are similar to std::integral_constant in that they
    provide value, value_type (and name).

SourceFiles
    exprTraits.C

Namespace
    Foam::expressions

Description
    A namespace for expression-related classes/traits etc.

Namespace
    Foam::expressions::Detail

Description
    A namespace for implementation details related to expressions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprTraits_H
#define Foam_expressions_exprTraits_H

// Regular field types
#include "label.H"
#include "scalar.H"
//TBD: #include "complex.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

//- An enumeration of known and expected expression value types.
//  Do not rely on the enumeration values for any direct coding.
//
//  \note NONE used when initializing types, whereas INVALID is used
//  for unsupported types (never as a stored type).
//  This avoids false positives when testing.
//
//  Except NONE and INVALID, the enumerations will mostly not be used
//  directly, but through exprTypeTraits :: value

enum class valueTypeCode : unsigned char
{
    NONE = 0,               //!< No type, or default initialized type
    INVALID,                //!< Invalid/unknown/error type

    // Rank 0 types
    type_bool,              //!< Type is 'bool'
    type_label,             //!< Type is 'label'
    type_scalar,            //!< Type is 'scalar'
//TBD: type_complex,           //!< Type is 'complex'

    // Rank 1 types
    type_vector,            //!< Type is 'vector'

    // Rank 2 types
    type_sphericalTensor,   //!< Type is 'sphericalTensor'
    type_symmTensor,        //!< Type is 'symmTensor'
    type_tensor             //!< Type is 'tensor'
};


// Global Functions

//- From string to valueTypeCode (if any)
valueTypeCode valueTypeCodeOf
(
    const word& dataTypeName,
    // Fallback for unknown
    const expressions::valueTypeCode deflt = expressions::valueTypeCode::INVALID
);


// Some implementation detail
namespace Detail
{

//- The number of components associated with given valueTypeCode
::Foam::direction nComponents(const expressions::valueTypeCode) noexcept;

//- The vector-space rank associated with given valueTypeCode
::Foam::direction rank(const expressions::valueTypeCode) noexcept;

} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions


/*---------------------------------------------------------------------------*\
                       Class exprTypeTraits Declaration
\*---------------------------------------------------------------------------*/

// Generic enumerated traits is INVALID (unknown)
template<class Type>
struct exprTypeTraits
{
    // The value type
    typedef Type value_type;
    // The type name (similar to pTraits typeName)
    static constexpr const char* const name = "";
    // The enumeration number associated with the type
    static constexpr
        ::Foam::expressions::valueTypeCode value =
        ::Foam::expressions::valueTypeCode::INVALID;
    // The rank of the type
    static constexpr ::Foam::direction rank = 0;
    // The number of components
    static constexpr ::Foam::direction nComponents = 0;
};


#undef  defineExprTypeTraits
#define defineExprTypeTraits(Type, Name, Rank, NumCmpts)                      \
    /*! \brief \c Name expression traits */                                   \
    template<>                                                                \
    struct exprTypeTraits<Type>                                               \
    {                                                                         \
        typedef Type value_type;                                              \
        static constexpr const char* const name = #Name ;                     \
        static constexpr                                                      \
            ::Foam::expressions::valueTypeCode value =                        \
            ::Foam::expressions::valueTypeCode::type_##Name ;                 \
        static constexpr ::Foam::direction rank = Rank ;                      \
        static constexpr ::Foam::direction nComponents = NumCmpts ;           \
    };


// Define with name to match pTraits::typeName, with rank/nComponents
defineExprTypeTraits(bool, bool, 0, 1);
defineExprTypeTraits(::Foam::label, label, 0, 1);
defineExprTypeTraits(::Foam::scalar, scalar, 0, 1);
//TBD: defineExprTypeTraits(::Foam::complex, complex, 0, 2);
defineExprTypeTraits(::Foam::vector, vector, 1, 3);
defineExprTypeTraits(::Foam::tensor, tensor, 2, 9);
defineExprTypeTraits(::Foam::sphericalTensor, sphericalTensor, 2, 1);
defineExprTypeTraits(::Foam::symmTensor, symmTensor, 2, 6);

#undef defineExprTypeTraits

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A word representation of a valueTypeCode.
//- Empty for expressions::valueTypeCode::INVALID
word name(const expressions::valueTypeCode typeCode);


//- A word representation of a valueTypeCode.
//- Empty for expressions::valueTypeCode::INVALID
template<>
struct nameOp<expressions::valueTypeCode>
{
    word operator()(const expressions::valueTypeCode typeCode) const
    {
        return Foam::name(typeCode);
    }
};

// No IOstream Operators for valueTypeCode at the moment (Nov 2021)

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
