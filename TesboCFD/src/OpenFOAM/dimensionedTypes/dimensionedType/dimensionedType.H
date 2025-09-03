/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::dimensioned

Description
    Generic dimensioned Type class

SourceFiles
    dimensionedType.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dimensionedType_H
#define Foam_dimensionedType_H

#include "word.H"
#include "direction.H"
#include "dimensionSet.H"
#include "VectorSpace.H"
#include "IOobjectOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class one;
class zero;
class dictionary;
class primitiveEntry;

template<class Type> class dimensioned;

template<class Type>
Istream& operator>>(Istream& is, dimensioned<Type>& dt);

/*---------------------------------------------------------------------------*\
                           Class dimensioned Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dimensioned
{
    // Private Data

        //- The variable name
        word name_;

        //- The dimension set
        dimensionSet dimensions_;

        //- The data value
        Type value_;


    // Private Member Functions

        //- Read helper.
        //  Requires a value, optional preceded with name and/or dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the name is present, it is used to rename.
        //  If dimensions are present, they are read.
        //  With checkDims = true, the dimensions read are verified
        //  against the current (expected) dimensions.
        void initialize(Istream& is, const bool checkDims);

        //- Find entry and assign to dimensioned Type
        //- FatalIOError if it is found and the number of tokens is incorrect,
        //- or it is mandatory and not found.
        //
        //  Requires a value, optional preceded with name and/or dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the name is present, it is used to rename.
        //  If dimensions are present, they are read.
        //  With checkDims = true, the dimensions read are verified
        //  against the current (expected) dimensions.
        bool readEntry
        (
            const word& key,
            const dictionary& dict,
            IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ,
            const bool checkDims = true, //!< verify dimensions read
            enum keyType::option matchOpt = keyType::REGEX
        );

public:

    //- The underlying data type
    typedef Type value_type;

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Constructors

        //- A dimensionless Zero, named "0"
        dimensioned();

        //- A dimensioned Zero, named "0"
        explicit dimensioned(const dimensionSet& dims);

        //- A dimensioned Zero, named "0"
        explicit dimensioned(const dimensionSet& dims, const Foam::zero);

        //- A dimensioned pTraits::one, named "1"
        explicit dimensioned(const dimensionSet& dims, const Foam::one);

        //- Implicit construct dimensionless from given value.
        dimensioned(const Type& val)
        :
            name_(::Foam::name(val)),
            dimensions_(dimless),
            value_(val)
        {}

        //- Construct dimensioned from given value
        dimensioned(const dimensionSet& dims, const Type& val);

        //- Construct from components (name, dimensions, value).
        dimensioned
        (
            const word& name,
            const dimensionSet& dims,
            const Type& val
        );

        //- Copy construct dimensioned Type with a new name
        dimensioned(const word& name, const dimensioned<Type>& dt);

        //- Construct from primitive entry with given name.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the optional name is found, it is used for renaming.
        //  If the optional dimensions are present, they are read and
        //  used without further verification.
        //  If no dimensions are found, the quantity is dimensionless.
        //  Fatal if not primitiveEntry or if number of tokens is incorrect.
        explicit dimensioned(const primitiveEntry& e);

        //- Construct from primitive entry with given name and dimensions.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the optional name is found, it is used for renaming.
        //  If the optional dimensions are present, they are read and
        //  verified against the expected dimensions.
        //  Fatal if not primitiveEntry or if number of tokens is incorrect.
        explicit dimensioned(const primitiveEntry& e, const dimensionSet& dims);

        //- Construct from dictionary lookup with a given name.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the optional name is found, it is used for renaming.
        //  If the optional dimensions are present, they are read and
        //  used without further verification.
        //  If no dimensions are found, the quantity is dimensionless.
        dimensioned(const word& name, const dictionary& dict);

        //- Construct from dictionary lookup with a given name and dimensions.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the optional name is found, it is used for renaming.
        //  If the optional dimensions are present, they are read and
        //  verified against the expected dimensions.
        dimensioned
        (
            const word& name,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct from dictionary lookup with a given name and dimensions.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        //  If the optional name is found, it is used for renaming.
        //  If the optional dimensions are present, they are read and
        //  verified against the expected dimensions.
        dimensioned
        (
            const word& name,
            const dimensionSet& dims,
            const dictionary& dict,
            const word& entryName      //!< dictionary lookup name
        );

        //- Construct from components (name, dimensions, value) with
        //- optional dictionary override.
        //  The entry may contain optional name and dimensions.
        //  \verbatim
        //  [name] [dims] value
        //  \endverbatim
        dimensioned
        (
            const word& name,
            const dimensionSet& dims,
            const Type& val,
            const dictionary& dict
        );


    // Static Member Functions

        //- Construct dimensioned from dictionary, with default value.
        //  FatalIOError if there are excess tokens.
        static dimensioned<Type> getOrDefault
        (
            const word& name,
            const dictionary& dict,
            const dimensionSet& dims = dimless,
            const Type& deflt = Type(Zero)
        );

        //- Construct dimensionless from dictionary, with default value.
        //  FatalIOError if it is found and there are excess tokens.
        static dimensioned<Type> getOrDefault
        (
            const word& name,
            const dictionary& dict,
            const Type& deflt = Type(Zero)
        );

        //- Construct dimensioned from dictionary, with default value.
        //  If the value is not found, it is added into the dictionary.
        //  FatalIOError if it is found and there are excess tokens.
        static dimensioned<Type> getOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const dimensionSet& dims = dimless,
            const Type& deflt = Type(Zero)
        );

        //- Construct dimensionless from dictionary, with default value.
        //  If the value is not found, it is added into the dictionary.
        //  FatalIOError if it is found and there are excess tokens.
        static dimensioned<Type> getOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const Type& deflt = Type(Zero)
        );


    // Member Functions

        //- Return const reference to name.
        const word& name() const noexcept { return name_; }

        //- Return non-const reference to name.
        word& name() noexcept { return name_; }

        //- Return const reference to dimensions.
        const dimensionSet& dimensions() const noexcept { return dimensions_; }

        //- Return non-const reference to dimensions.
        dimensionSet& dimensions() noexcept { return dimensions_; }

        //- Return const reference to value.
        const Type& value() const noexcept { return value_; }

        //- Return non-const reference to value.
        Type& value() noexcept { return value_; }

        //- Return a component as a dimensioned<cmptType>
        dimensioned<cmptType> component(const direction d) const;

        //- Return a component with a dimensioned<cmptType>
        void replace(const direction d, const dimensioned<cmptType>& dc);

        //- Return transpose.
        dimensioned<Type> T() const;

        //- Update the value of dimensioned\<Type\>,
        //- lookup in dictionary with the name().
        bool read(const dictionary& dict);

        //- Update the value of dimensioned\<Type\> if found in the dictionary,
        //- lookup in dictionary with the name().
        bool readIfPresent(const dictionary& dict);

        //- Update the value of dimensioned\<Type\>,
        //- using an alternative entry name
        bool read(const word& entryName, const dictionary& dict);

        //- Update the value of dimensioned\<Type\> if found in the dictionary,
        //- using an alternative entry name
        bool readIfPresent(const word& entryName, const dictionary& dict);


    // IO

        //- Read (name, dimensions, value) from stream,
        //- using units from system table.
        //  Optionally skip reading the name
        Istream& read(Istream& is, const bool readName = true);

        //- Read (name, dimensions, value) from stream,
        //- using units from dictionary
        Istream& read(Istream& is, const dictionary& readSet);

        //- Read (name, dimensions, value) from stream,
        //- using units from table
        Istream& read(Istream& is, const HashTable<dimensionedScalar>& readSet);

        //- Write as a dictionary entry with keyword.
        //  The name is not written when it is identical to keyword.
        //  The dimensions are always written.
        void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- Return a component as a dimensioned<cmptType>
        dimensioned<cmptType> operator[](const direction d) const;

        void operator+=(const dimensioned<Type>& dt);
        void operator-=(const dimensioned<Type>& dt);
        void operator*=(const scalar s);
        void operator/=(const scalar s);


    // IOstream Operators

        //- Read from stream. The name and dimensions are optional.
        //  If the optional dimensions are present,
        //  they are used without further verification.
        friend Istream& operator>> <Type>
        (
            Istream& is,
            dimensioned<Type>& dt
        );


    // Housekeeping

        //- Deprecated(2018-11) Construct from Istream
        //- (expects name, dimensions, value)
        //  \deprecated(2018-11) - should generally use construct from
        //      dictionary or primitiveEntry instead
        //      (additional checks on the input stream).
        FOAM_DEPRECATED_FOR(2018-11, "construct from dictionary or entry")
        explicit dimensioned(Istream& is);

        //- Deprecated(2018-11) Construct from Istream with given name
        //- (expects dimensions, value)
        //  \deprecated(2018-11) - should generally use construct from
        //      dictionary or primitiveEntry instead
        //      (additional checks on the input stream).
        FOAM_DEPRECATED_FOR(2018-11, "construct from dictionary or entry")
        dimensioned(const word& name, Istream& is);

        //- Deprecated(2018-11) Construct from Istream with given name
        //- and expected dimensions.
        //  Expects value, but supports optional name and dimensions.
        //  If the optional dimensions are present, they are read and
        //  verified against the expected dimensions.
        //  \deprecated(2018-11) - should generally use construct from
        //      dictionary or primitiveEntry instead
        //      (additional checks on the input stream).
        FOAM_DEPRECATED_FOR(2018-11, "construct from dictionary or entry")
        dimensioned(const word& name, const dimensionSet& dims, Istream& is);

        //- Same as getOrDefault()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrDefault()")
        static dimensioned<Type> lookupOrDefault
        (
            const word& name,
            const dictionary& dict,
            const dimensionSet& dims = dimless,
            const Type& deflt = Type(Zero)
        )
        {
            return getOrDefault(name, dict, dims, deflt);
        }

        //- Same as getOrDefault()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrDefault()")
        static dimensioned<Type> lookupOrDefault
        (
            const word& name,
            const dictionary& dict,
            const Type& deflt = Type(Zero)
        )
        {
            return getOrDefault(name, dict, deflt);
        }

        //- Same as getOrAddToDict()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrAddToDict()")
        static dimensioned<Type> lookupOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const dimensionSet& dims = dimless,
            const Type& deflt = Type(Zero)
        )
        {
            return getOrAddToDict(name, dict, dims, deflt);
        }

        //- Same as getOrAddToDict()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrAddToDict()")
        static dimensioned<Type> lookupOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const Type& deflt = Type(Zero)
        )
        {
            return getOrAddToDict(name, dict, deflt);
        }
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Output operator
template<class Type>
Ostream& operator<<(Ostream& os, const dimensioned<Type>& dt);

template<class Type, direction r>
dimensioned<typename powProduct<Type, r>::type>
pow
(
    const dimensioned<Type>&,
    typename powProduct<Type, r>::type
  = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type>
dimensioned<typename outerProduct<Type, Type>::type>
sqr(const dimensioned<Type>&);

template<class Type>
dimensioned<typename typeOfMag<Type>::type>
magSqr(const dimensioned<Type>& dt);

template<class Type>
dimensioned<typename typeOfMag<Type>::type>
mag(const dimensioned<Type>& dt);

template<class Type>
dimensioned<Type> cmptMultiply
(
    const dimensioned<Type>&,
    const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> cmptDivide
(
    const dimensioned<Type>&,
    const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> max(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> min(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> lerp
(
    const dimensioned<Type>& a,
    const dimensioned<Type>& b,
    const scalar t
);

template<class Type>
bool operator<(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
bool operator>(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator+(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator-(const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator-(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator*
(
    const dimensioned<scalar>&,
    const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> operator/
(
    const dimensioned<Type>&,
    const dimensioned<scalar>&
);


#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
dimensioned<typename product<Type1, Type2>::type>                              \
operator op(const dimensioned<Type1>&, const dimensioned<Type2>&);             \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
dimensioned<typename product<Type, Form>::type>                                \
operator op                                                                    \
(                                                                              \
    const dimensioned<Type>&,                                                  \
    const VectorSpace<Form,Cmpt,nCmpt>&                                        \
);                                                                             \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
dimensioned<typename product<Form, Type>::type>                                \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>&,                                       \
    const dimensioned<Type>&                                                   \
);

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dimensionedType.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
