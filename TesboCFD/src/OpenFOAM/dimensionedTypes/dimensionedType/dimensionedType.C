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

\*---------------------------------------------------------------------------*/

#include "dimensionedType.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::dimensioned<Type>::initialize(Istream& is, const bool checkDims)
{
    token nextToken(is);
    is.putBack(nextToken);

    // Optional name found - use it
    if (nextToken.isWord())
    {
        is >> name_;
        is >> nextToken;
        is.putBack(nextToken);
    }

    scalar mult{1};

    if (nextToken == token::BEGIN_SQR)
    {
        // Optional dimensions found - use them
        const dimensionSet curr(dimensions_);
        dimensions_.read(is, mult);

        if (checkDims && curr != dimensions_)
        {
            FatalIOErrorInFunction(is)
                << "The dimensions " << dimensions_
                << " provided do not match the expected dimensions "
                << curr << endl
                << abort(FatalIOError);
        }
    }

    // Read value
    is >> value_;
    value_ *= mult;
}


template<class Type>
bool Foam::dimensioned<Type>::readEntry
(
    const word& key,
    const dictionary& dict,
    IOobjectOption::readOption readOpt,
    const bool checkDims,
    enum keyType::option matchOpt
)
{
    if (readOpt == IOobjectOption::NO_READ)
    {
        return false;
    }

    // Largely identical to dictionary::readEntry(),
    // but with optional handling of checkDims

    const entry* eptr = dict.findEntry(key, matchOpt);

    if (eptr)
    {
        ITstream& is = eptr->stream();

        initialize(is, checkDims);

        dict.checkITstream(is, key);

        return true;
    }
    else if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << key << "' not found in dictionary "
            << dict.name()
            << exit(FatalIOError);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<Type>::dimensioned()
:
    name_("0"),
    dimensions_(),
    value_(Zero)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned(const dimensionSet& dims)
:
    name_("0"),
    dimensions_(dims),
    value_(Zero)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const dimensionSet& dims,
    const Foam::zero
)
:
    name_("0"),
    dimensions_(dims),
    value_(Zero)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const dimensionSet& dims,
    const Foam::one
)
:
    name_("1"),
    dimensions_(dims),
    value_(pTraits<Type>::one)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const dimensionSet& dims,
    const Type& val
)
:
    name_(::Foam::name(val)),
    dimensions_(dims),
    value_(val)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dims,
    const Type& val
)
:
    name_(name),
    dimensions_(dims),
    value_(val)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensioned<Type>& dt
)
:
    name_(name),
    dimensions_(dt.dimensions_),
    value_(dt.value_)
{}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const primitiveEntry& e
)
:
    name_(e.name()),
    dimensions_(),
    value_(Zero)
{
    ITstream& is = e.stream();

    // checkDims = false
    initialize(is, false);

    e.checkITstream(is);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const primitiveEntry& e,
    const dimensionSet& dims
)
:
    name_(e.name()),
    dimensions_(dims),
    value_(Zero)
{
    ITstream& is = e.stream();

    // checkDims = true
    initialize(is, true);

    e.checkITstream(is);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dimensions_(),
    value_(Zero)
{
    // checkDims = false
    readEntry(name, dict, IOobjectOption::MUST_READ, false);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    name_(name),
    dimensions_(dims),
    value_(Zero)
{
    // checkDims = true
    readEntry(name, dict, IOobjectOption::MUST_READ);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict,
    const word& entryName
)
:
    name_(name),
    dimensions_(dims),
    value_(Zero)
{
    // checkDims = true
    readEntry(entryName, dict, IOobjectOption::MUST_READ);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dims,
    const Type& val,
    const dictionary& dict
)
:
    name_(name),
    dimensions_(dims),
    value_(val)
{
    // checkDims = true
    readEntry(name, dict, IOobjectOption::READ_IF_PRESENT);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned(Istream& is)
:
    dimensions_()
{
    read(is);
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    Istream& is
)
:
    name_(name),
    dimensions_()
{
    read(is, false);  // Don't read name. Read dimensionSet + multiplier only.
}


template<class Type>
Foam::dimensioned<Type>::dimensioned
(
    const word& name,
    const dimensionSet& dims,
    Istream& is
)
:
    name_(name),
    dimensions_(dims),
    value_(Zero)
{
    // checkDims = true
    initialize(is, true);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::getOrDefault
(
    const word& name,
    const dictionary& dict,
    const dimensionSet& dims,
    const Type& deflt
)
{
    // checkDims = true
    return dimensioned<Type>(name, dims, deflt, dict);
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::getOrDefault
(
    const word& name,
    const dictionary& dict,
    const Type& deflt
)
{
    return dimensioned<Type>(name, dimless, deflt, dict);
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::getOrAddToDict
(
    const word& name,
    dictionary& dict,
    const dimensionSet& dims,
    const Type& deflt
)
{
    if (dict.found(name))
    {
        return dimensioned<Type>(name, dims, dict);
    }

    (void) dict.add(name, deflt);
    return dimensioned<Type>(name, dims, deflt);
}


template<class Type>
Foam::dimensioned<Type> Foam::dimensioned<Type>::getOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Type& deflt
)
{
    return getOrAddToDict(name, dict, dimless, deflt);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<typename Foam::dimensioned<Type>::cmptType>
Foam::dimensioned<Type>::component
(
    const direction d
) const
{
    return dimensioned<cmptType>
    (
        name_ + ".component(" + Foam::name(d) + ')',
        dimensions_,
        value_.component(d)
    );
}


template<class Type>
void Foam::dimensioned<Type>::replace
(
    const direction d,
    const dimensioned<typename dimensioned<Type>::cmptType>& dc
)
{
    dimensions_ = dc.dimensions();
    value_.replace(d, dc.value());
}


template<class Type>
bool Foam::dimensioned<Type>::read(const dictionary& dict)
{
    return read(name_, dict);
}


template<class Type>
bool Foam::dimensioned<Type>::readIfPresent(const dictionary& dict)
{
    return readIfPresent(name_, dict);
}


template<class Type>
bool Foam::dimensioned<Type>::read
(
    const word& entryName,
    const dictionary& dict
)
{
    // checkDims = true
    return readEntry(entryName, dict, IOobjectOption::MUST_READ);
}


template<class Type>
bool Foam::dimensioned<Type>::readIfPresent
(
    const word& entryName,
    const dictionary& dict
)
{
    // checkDims = true
    return readEntry(entryName, dict, IOobjectOption::READ_IF_PRESENT);
}


template<class Type>
Foam::Istream& Foam::dimensioned<Type>::read(Istream& is, const bool readName)
{
    if (readName)
    {
        // Read name
        is >> name_;
    }

    // Read dimensionSet + multiplier
    scalar mult{1};
    dimensions_.read(is, mult);

    // Read value
    is >> value_;
    value_ *= mult;

    is.check(FUNCTION_NAME);
    return is;
}


template<class Type>
Foam::Istream&
Foam::dimensioned<Type>::read(Istream& is, const dictionary& readSet)
{
    // Read name
    is >> name_;

    // Read dimensionSet + multiplier
    scalar mult{1};
    dimensions_.read(is, mult, readSet);

    // Read value
    is >> value_;
    value_ *= mult;

    is.check(FUNCTION_NAME);
    return is;
}


template<class Type>
Foam::Istream& Foam::dimensioned<Type>::read
(
    Istream& is,
    const HashTable<dimensionedScalar>& readSet
)
{
    // Read name
    is >> name_;

    // Read dimensionSet + multiplier
    scalar mult{1};
    dimensions_.read(is, mult, readSet);

    // Read value
    is >> value_;
    value_ *= mult;

    is.check(FUNCTION_NAME);
    return is;
}


template<class Type>
void Foam::dimensioned<Type>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);

        if (keyword != name_)
        {
            // The name, only if different from keyword
            os << name_ << token::SPACE;
        }
    }
    else
    {
        // Use the name (no keyword provided)
        os.writeKeyword(name_);
    }

    // The dimensions
    scalar mult{1};
    dimensions_.write(os, mult);

    // The value
    os << token::SPACE << value_/mult;
    os.endEntry();

    os.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::dimensioned<typename Foam::dimensioned<Type>::cmptType>
Foam::dimensioned<Type>::operator[]
(
    const direction d
) const
{
    return component(d);
}


template<class Type>
void Foam::dimensioned<Type>::operator+=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ += dt.dimensions_;
    value_ += dt.value_;
}


template<class Type>
void Foam::dimensioned<Type>::operator-=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ -= dt.dimensions_;
    value_ -= dt.value_;
}


template<class Type>
void Foam::dimensioned<Type>::operator*=
(
    const scalar s
)
{
    value_ *= s;
}


template<class Type>
void Foam::dimensioned<Type>::operator/=
(
    const scalar s
)
{
    value_ /= s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type, Foam::direction r>
Foam::dimensioned<typename Foam::powProduct<Type, r>::type>
Foam::pow(const dimensioned<Type>& dt, typename powProduct<Type, r>::type)
{
    return dimensioned<typename powProduct<Type, r>::type>
    (
        "pow(" + dt.name() + ',' + name(r) + ')',
        pow(dt.dimensions(), r),
        pow(dt.value(), 2)
    );
}


template<class Type>
Foam::dimensioned<typename Foam::outerProduct<Type, Type>::type>
Foam::sqr(const dimensioned<Type>& dt)
{
    return dimensioned<typename outerProduct<Type, Type>::type>
    (
        "sqr(" + dt.name() + ')',
        sqr(dt.dimensions()),
        sqr(dt.value())
    );
}

template<class Type>
Foam::dimensioned<typename Foam::typeOfMag<Type>::type>
Foam::magSqr(const dimensioned<Type>& dt)
{
    typedef typename typeOfMag<Type>::type magType;

    return dimensioned<magType>
    (
        "magSqr(" + dt.name() + ')',
        magSqr(dt.dimensions()),
        magSqr(dt.value())
    );
}

template<class Type>
Foam::dimensioned<typename Foam::typeOfMag<Type>::type>
Foam::mag(const dimensioned<Type>& dt)
{
    typedef typename typeOfMag<Type>::type magType;

    return dimensioned<magType>
    (
        "mag(" + dt.name() + ')',
        dt.dimensions(),
        mag(dt.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::cmptMultiply
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        "cmptMultiply(" + dt1.name() + ',' + dt2.name() + ')',
        cmptMultiply(dt1.dimensions(), dt2.dimensions()),
        cmptMultiply(dt1.value(), dt2.value())
    );
}

template<class Type>
Foam::dimensioned<Type> Foam::cmptDivide
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        "cmptDivide(" + dt1.name() + ',' + dt2.name() + ')',
        cmptDivide(dt1.dimensions(), dt2.dimensions()),
        cmptDivide(dt1.value(), dt2.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::max
(
    const dimensioned<Type>& a,
    const dimensioned<Type>& b
)
{
    return dimensioned<Type>
    (
        "max(" + a.name() + ',' + b.name() + ')',
        max(a.dimensions(), b.dimensions()),
        max(a.value(), b.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::min
(
    const dimensioned<Type>& a,
    const dimensioned<Type>& b
)
{
    return dimensioned<Type>
    (
        "min(" + a.name() + ',' + b.name() + ')',
        min(a.dimensions(), b.dimensions()),
        min(a.value(), b.value())
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::lerp
(
    const dimensioned<Type>& a,
    const dimensioned<Type>& b,
    const scalar t
)
{
    return dimensioned<Type>
    (
        "lerp(" + a.name() + ',' + b.name() + ',' + ::Foam::name(t) + ')',
        lerp(a.dimensions(), b.dimensions()),
        lerp(a.value(), b.value(), t)
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, dimensioned<Type>& dt)
{
    dt.initialize(is, false);  // no checkDims
    is.check(FUNCTION_NAME);
    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const dimensioned<Type>& dt)
{
    // The name
    os << dt.name() << token::SPACE;

    // The dimensions
    scalar mult{1};
    dt.dimensions().write(os, mult);

    // The value
    os << token::SPACE << dt.value()/mult;

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
bool Foam::operator<
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dt1.value() < dt2.value();
}


template<class Type>
bool Foam::operator>
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dt2.value() < dt1.value();
}


template<class Type>
Foam::dimensioned<Type> Foam::operator+
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        '(' + dt1.name() + '+' + dt2.name() + ')',
        dt1.dimensions() + dt2.dimensions(),
        dt1.value() + dt2.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator-(const dimensioned<Type>& dt)
{
    return dimensioned<Type>
    (
        '-' + dt.name(),
        dt.dimensions(),
        -dt.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator-
(
    const dimensioned<Type>& dt1,
    const dimensioned<Type>& dt2
)
{
    return dimensioned<Type>
    (
        '(' + dt1.name() + '-' + dt2.name() + ')',
        dt1.dimensions() - dt2.dimensions(),
        dt1.value() - dt2.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const dimensioned<Type>& dt
)
{
    return dimensioned<Type>
    (
        '(' + ds.name() + '*' + dt.name() + ')',
        ds.dimensions() * dt.dimensions(),
        ds.value() * dt.value()
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::operator/
(
    const dimensioned<Type>& dt,
    const dimensioned<scalar>& ds
)
{
    return dimensioned<Type>
    (
        '(' + dt.name() + '|' + ds.name() + ')',
        dt.dimensions() / ds.dimensions(),
        dt.value() / ds.value()
    );
}


#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
Foam::dimensioned<typename Foam::product<Type1, Type2>::type>                  \
Foam::operator op                                                              \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Type1, Type2>::type>                   \
    (                                                                          \
        '(' + dt1.name() + #op + dt2.name() + ')',                             \
        dt1.dimensions() op dt2.dimensions(),                                  \
        dt1.value() op dt2.value()                                             \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, Foam::direction nCmpt>            \
Foam::dimensioned<typename Foam::product<Type, Form>::type>                    \
Foam::operator op                                                              \
(                                                                              \
    const dimensioned<Type>& dt1,                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& t2                                     \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Type, Form>::type>                     \
    (                                                                          \
        '(' + dt1.name() + #op + name(t2) + ')',                               \
        dt1.dimensions(),                                                      \
        dt1.value() op static_cast<const Form&>(t2)                            \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, Foam::direction nCmpt>            \
Foam::dimensioned<typename Foam::product<Form, Type>::type>                    \
Foam::operator op                                                              \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& t1,                                    \
    const dimensioned<Type>& dt2                                               \
)                                                                              \
{                                                                              \
    return dimensioned<typename product<Form, Type>::type>                     \
    (                                                                          \
        '(' + name(t1) + #op + dt2.name() + ')',                               \
        dt2.dimensions(),                                                      \
        static_cast<const Form&>(t1) op dt2.value()                            \
    );                                                                         \
}


PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// ************************************************************************* //
