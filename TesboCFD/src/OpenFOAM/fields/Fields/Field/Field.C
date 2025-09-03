/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "FieldMapper.H"
#include "FieldM.H"
#include "dictionary.H"
#include "contiguous.H"
#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const labelUList& mapAddressing
)
:
    List<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
:
    List<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing);
}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
:
    List<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
:
    List<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const FieldMapper& mapper,
    const bool applyFlip
)
:
    List<Type>(mapper.size())
{
    map(mapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const FieldMapper& mapper,
    const Type& defaultValue,
    const bool applyFlip
)
:
    List<Type>(mapper.size(), defaultValue)
{
    map(mapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const FieldMapper& mapper,
    const UList<Type>& defaultValues,
    const bool applyFlip
)
:
    List<Type>(defaultValues)
{
    map(mapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const FieldMapper& mapper,
    const bool applyFlip
)
:
    List<Type>(mapper.size())
{
    map(tmapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const FieldMapper& mapper,
    const Type& defaultValue,
    const bool applyFlip
)
:
    List<Type>(mapper.size(), defaultValue)
{
    map(tmapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const FieldMapper& mapper,
    const UList<Type>& defaultValues,
    const bool applyFlip
)
:
    List<Type>(defaultValues)
{
    map(tmapF, mapper, applyFlip);
}


template<class Type>
Foam::Field<Type>::Field(const entry& e, const label len)
{
    Field<Type>::assign(e, len);
}


template<class Type>
Foam::Field<Type>::Field
(
    const word& key,
    const dictionary& dict,
    const label len,
    IOobjectOption::readOption readOpt
)
{
    if (!Field<Type>::assign(key, dict, len, readOpt))
    {
        if (IOobjectOption::isReadOptional(readOpt))
        {
            // Lazy read: init with zero value
            if (len > 0) this->resize(len, Zero);
        }
        else
        {
            // No read: set length only
            if (len > 0) this->resize(len);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Field<Type>::assign(const entry& e, const label len)
{
    if (len)
    {
        ITstream& is = e.stream();

        // Read first token
        token firstToken(is);

        if (firstToken.isWord("uniform"))
        {
            // Resize to expected length (or -1 : retain current length)
            if (len >= 0)
            {
                this->resize_nocopy(len);
            }
            operator=(pTraits<Type>(is));
        }
        else if (firstToken.isWord("nonuniform"))
        {
            is >> static_cast<List<Type>&>(*this);
            const label lenRead = this->size();

            // Check lengths
            if (len >= 0 && len != lenRead)
            {
                if (len < lenRead && allowConstructFromLargerSize)
                {
                    // Truncate the data
                    this->resize(len);

                    #ifdef FULLDEBUG
                    IOWarningInFunction(is)
                        << "Sizes do not match. Truncating " << lenRead
                        << " entries to " << len << endl;
                    #endif
                }
                else
                {
                    FatalIOErrorInFunction(is)
                        << "Size " << lenRead
                        << " is not equal to the expected length " << len
                        << exit(FatalIOError);
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Expected keyword 'uniform' or 'nonuniform', found "
                << firstToken.info() << nl
                << exit(FatalIOError);
        }
    }
}


template<class Type>
bool Foam::Field<Type>::assign
(
    const word& key,
    const dictionary& dict,
    const label len,
    IOobjectOption::readOption readOpt
)
{
    if (!len)
    {
        return true;
    }
    else if (readOpt != IOobjectOption::NO_READ)
    {
        const entry* eptr = dict.findEntry(key, keyType::LITERAL);

        if (eptr)
        {
            Field<Type>::assign(*eptr, len);
            return true;
        }

        // Missing (mandatory or optional)
        if (IOobjectOption::isReadRequired(readOpt))
        {
            FatalIOErrorInFunction(dict)
                << "Required entry '" << key << "' missing in dictionary "
                << dict.relativeName() << nl
                << exit(FatalIOError);
        }
    }

    return false;
}


template<class Type>
void Foam::Field<Type>::map
(
    const UList<Type>& mapF,
    const labelUList& mapAddressing
)
{
    Field<Type>& f = *this;

    if (f.size() != mapAddressing.size())
    {
        f.resize(mapAddressing.size());
    }

    if (mapF.size() > 0)
    {
        forAll(f, i)
        {
            const label mapI = mapAddressing[i];

            if (mapI >= 0)
            {
                f[i] = mapF[mapI];
            }
        }
    }
}


template<class Type>
void Foam::Field<Type>::map
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
{
    map(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::map
(
    const UList<Type>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
{
    Field<Type>& f = *this;

    if (f.size() != mapAddressing.size())
    {
        f.resize(mapAddressing.size());
    }

    if (mapWeights.size() != mapAddressing.size())
    {
        FatalErrorInFunction
            << mapWeights.size() << " map size: " << mapAddressing.size()
            << abort(FatalError);
    }

    forAll(f, i)
    {
        const labelList&  localAddrs   = mapAddressing[i];
        const scalarList& localWeights = mapWeights[i];

        f[i] = Zero;

        forAll(localAddrs, j)
        {
            f[i] += localWeights[j]*mapF[localAddrs[j]];
        }
    }
}


template<class Type>
void Foam::Field<Type>::map
(
    const tmp<Field<Type>>& tmapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
{
    map(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::map
(
    const UList<Type>& mapF,
    const FieldMapper& mapper,
    const bool applyFlip
)
{
    if (mapper.distributed())
    {
        // Fetch remote parts of mapF
        const mapDistributeBase& distMap = mapper.distributeMap();
        Field<Type> newMapF(mapF);

        if (applyFlip)
        {
            distMap.distribute(newMapF);
        }
        else
        {
            distMap.distribute(newMapF, identityOp());
        }

        if (mapper.direct() && notNull(mapper.directAddressing()))
        {
            map(newMapF, mapper.directAddressing());
        }
        else if (!mapper.direct())
        {
            map(newMapF, mapper.addressing(), mapper.weights());
        }
        else if (mapper.direct() && isNull(mapper.directAddressing()))
        {
            // Special case, no local mapper. Assume ordering already correct
            // from distribution. Note: this behaviour is different compared
            // to local mapper.
            this->transfer(newMapF);
            this->resize(mapper.size());
        }
    }
    else
    {
        if
        (
            mapper.direct()
         && notNull(mapper.directAddressing())
         && mapper.directAddressing().size()
        )
        {
            map(mapF, mapper.directAddressing());
        }
        else if (!mapper.direct() && mapper.addressing().size())
        {
            map(mapF, mapper.addressing(), mapper.weights());
        }
    }
}


template<class Type>
void Foam::Field<Type>::map
(
    const tmp<Field<Type>>& tmapF,
    const FieldMapper& mapper,
    const bool applyFlip
)
{
    map(tmapF(), mapper, applyFlip);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::autoMap
(
    const FieldMapper& mapper,
    const bool applyFlip
)
{
    if (mapper.distributed())
    {
        // Fetch remote parts of *this
        const mapDistributeBase& distMap = mapper.distributeMap();
        Field<Type> fCpy(*this);

        if (applyFlip)
        {
            distMap.distribute(fCpy);
        }
        else
        {
            distMap.distribute(fCpy, identityOp());
        }

        if
        (
            (mapper.direct()
         && notNull(mapper.directAddressing()))
         || !mapper.direct()
        )
        {
            this->map(fCpy, mapper);
        }
        else if (mapper.direct() && isNull(mapper.directAddressing()))
        {
            // Special case, no local mapper. Assume ordering already correct
            // from distribution. Note: this behaviour is different compared
            // to local mapper.
            this->transfer(fCpy);
            this->resize(mapper.size());
        }
    }
    else
    {
        if
        (
            (
                mapper.direct()
             && notNull(mapper.directAddressing())
             && mapper.directAddressing().size()
            )
         || (!mapper.direct() && mapper.addressing().size())
        )
        {
            Field<Type> fCpy(*this);
            map(fCpy, mapper);
        }
        else
        {
            this->resize(mapper.size());
        }
    }
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const UList<Type>& mapF,
    const labelUList& mapAddressing
)
{
    Field<Type>& f = *this;

    forAll(mapF, i)
    {
        label mapI = mapAddressing[i];

        if (mapI >= 0)
        {
            f[mapI] = mapF[i];
        }
    }
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
{
    rmap(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const UList<Type>& mapF,
    const labelUList& mapAddressing,
    const UList<scalar>& mapWeights
)
{
    Field<Type>& f = *this;

    f = Zero;

    forAll(mapF, i)
    {
        f[mapAddressing[i]] += mapF[i]*mapWeights[i];
    }
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing,
    const UList<scalar>& mapWeights
)
{
    rmap(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::negate()
{
    TFOR_ALL_F_OP_OP_F(Type, *this, =, -, Type, *this)
}


// A no-op except for vector specialization
template<class Type>
void Foam::Field<Type>::normalise()
{}


template<class Type>
Foam::tmp<Foam::Field<typename Foam::Field<Type>::cmptType>>
Foam::Field<Type>::component
(
    const direction d
) const
{
    auto tres = tmp<Field<cmptType>>::New(this->size());
    ::Foam::component(tres.ref(), *this, d);
    return tres;
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const UList<cmptType>& sf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, *this, ., replace, const direction, d,
        cmptType, sf)
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const tmp<Field<cmptType>>& tsf
)
{
    replace(d, tsf());
    tsf.clear();
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const cmptType& c
)
{
    TFOR_ALL_F_OP_FUNC_S_S(Type, *this, ., replace, const direction, d,
        cmptType, c)
}


template<class Type>
void Foam::Field<Type>::clamp_min(const Type& lower)
{
    // Use free function max() [sic] to impose component-wise clamp_min
    // std::for_each
    for (auto& val : *this)
    {
        val = max(val, lower);
    }
}

template<class Type>
void Foam::Field<Type>::clamp_max(const Type& upper)
{
    // Use free function min() [sic] to impose component-wise clamp_max
    // std::for_each
    for (auto& val : *this)
    {
        val = min(val, upper);
    }
}


template<class Type>
void Foam::Field<Type>::clamp_range(const Type& lower, const Type& upper)
{
    // Note: no checks for bad/invalid clamping ranges

    // Use free functions min(), max() to impose component-wise clamping
    // std::for_each
    for (auto& val : *this)
    {
        val = min(max(val, lower), upper);
    }
}

template<class Type>
void Foam::Field<Type>::clamp_range(const MinMax<Type>& range)
{
    Field<Type>::clamp_range(range.min(), range.max());
}


template<class Type>
template<class VSForm>
VSForm Foam::Field<Type>::block(const label start) const
{
    VSForm vs;
    for (direction i=0; i<VSForm::nComponents; i++)
    {
        vs[i] = this->operator[](start + i);
    }
    return vs;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Field<Type>::T() const
{
    auto tres = tmp<Field<Type>>::New(this->size());
    ::Foam::T(tres.ref(), *this);
    return tres;
}


template<class Type>
void Foam::Field<Type>::writeEntry(const word& keyword, Ostream& os) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }

    // The contents are 'uniform' if the list is non-empty
    // and all entries have identical values.

    if (is_contiguous<Type>::value && List<Type>::uniform())
    {
        os << word("uniform") << token::SPACE << List<Type>::front();
    }
    else
    {
        os << word("nonuniform") << token::SPACE;
        List<Type>::writeEntry(os);
    }

    os.endEntry();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Field<Type>::operator=(const Field<Type>& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    List<Type>::operator=(rhs);
}


template<class Type>
void Foam::Field<Type>::operator=(const tmp<Field>& rhs)
{
    if (this == &(rhs()))
    {
        return;  // Self-assignment is a no-op
    }

    List<Type>::operator=(rhs());
}


template<class Type>
template<class Form, class Cmpt, Foam::direction nCmpt>
void Foam::Field<Type>::operator=(const VectorSpace<Form,Cmpt,nCmpt>& vs)
{
    TFOR_ALL_F_OP_S(Type, *this, =, VSType, vs)
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const UList<TYPE>& f)                      \
{                                                                              \
    TFOR_ALL_F_OP_F(Type, *this, op, TYPE, f)                                  \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const tmp<Field<TYPE>>& tf)                \
{                                                                              \
    operator op(tf());                                                         \
    tf.clear();                                                                \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const TYPE& t)                             \
{                                                                              \
    TFOR_ALL_F_OP_S(Type, *this, op, TYPE, t)                                  \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const Field<Type>& f)
{
    os  << static_cast<const List<Type>&>(f);
    return os;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const tmp<Field<Type>>& tf)
{
    os  << tf();
    tf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFunctions.C"

// ************************************************************************* //
