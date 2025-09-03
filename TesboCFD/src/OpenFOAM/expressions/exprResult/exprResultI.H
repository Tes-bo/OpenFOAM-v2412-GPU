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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
inline bool Foam::expressions::exprResult::deleteChecked()
{
    const bool ok = isType<Type>();

    if (ok && fieldPtr_ != nullptr)
    {
        delete static_cast<Field<Type>*>(fieldPtr_);
        fieldPtr_ = nullptr;
        size_ = 0;
    }

    return ok;
}


template<class Type>
inline bool Foam::expressions::exprResult::readChecked
(
    const entry& e,
    const label len,
    const bool singleValueOnly
)
{
    const bool ok = isType<Type>();

    if (ok)
    {
        destroy();

        if (singleValueOnly)
        {
            const Type val(e.get<Type>());

            size_ = len;
            fieldPtr_ = new Field<Type>(size_, val);

            value_.set(val);
        }
        else
        {
            size_ = len;
            fieldPtr_ = new Field<Type>(e, size_);

            value_.clear();
        }
    }

    return ok;
}


template<class Type>
bool Foam::expressions::exprResult::getUniformChecked
(
    exprResult& result,
    const label size,
    const bool noWarn,
    const bool parRun
) const
{
    if (!isType<Type>())
    {
        return false;
    }

    result.clear();

    const Field<Type>& fld = *static_cast<const Field<Type>*>(fieldPtr_);

    const Type avg = (parRun ? gAverage(fld) : average(fld));

    if (!noWarn)
    {
        const MinMax<Type> limits = (parRun ? gMinMax(fld) : minMax(fld));

        if (limits.mag() > SMALL)
        {
            WarningInFunction
                << "Different min/max values: " << limits
                << " Using the average " << avg << nl;
        }
    }

    result.setResult(avg, size);

    return true;
}


template<class Type>
bool Foam::expressions::exprResult::plusEqChecked
(
    const exprResult& b
)
{
    const bool ok = isType<Type>();

    if (ok)
    {
          *static_cast<Field<Type>*>(fieldPtr_)
       += *static_cast<const Field<Type>*>(b.fieldPtr_);
    }

    return ok;
}


template<class Type>
bool Foam::expressions::exprResult::multiplyEqChecked
(
    const scalar& b
)
{
    const bool ok = isType<Type>();

    if (ok)
    {
        *static_cast<Field<Type>*>(fieldPtr_) *= b;
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::expressions::exprResult::exprResult(const Field<Type>& fld)
:
    exprResult()
{
    DebugInFunction << nl;

    setResult(fld);
}


template<class Type>
Foam::expressions::exprResult::exprResult(Field<Type>&& fld)
:
    exprResult()
{
    DebugInFunction << nl;

    setResult(std::move(fld));
}


template<class Type>
Foam::expressions::exprResult::exprResult(const dimensioned<Type>& dt)
:
    exprResult()
{
    DebugInFunction << nl;
    setSingleValue(dt.value());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::exprResult::hasValue() const
{
    return (!valType_.empty() && fieldPtr_ != nullptr);
}


inline const Foam::word&
Foam::expressions::exprResult::valueType() const noexcept
{
    return valType_;
}


inline bool Foam::expressions::exprResult::isPointData
(
    const bool wantPointData
) const
{
    return isPointData_ == wantPointData;
}


inline bool Foam::expressions::exprResult::isUniform() const
{
    return value_.good();
}


template<class Type>
inline bool Foam::expressions::exprResult::isType() const
{
    return valType_ == pTraits<Type>::typeName;
}


template<class Type>
inline Type Foam::expressions::exprResult::getValue() const
{
    if (value_.good() && isType<Type>())
    {
        return value_.get<Type>();
    }

    // Fallthrough
    return Zero;
}


template<class Type>
inline const Foam::Field<Type>* Foam::expressions::exprResult::getField() const
{
    if (isType<Type>())
    {
        return static_cast<const Field<Type>*>(fieldPtr_);
    }

    // Fallthrough
    return nullptr;
}


inline bool Foam::expressions::exprResult::is_bool() const
{
    return valType_ == pTraits<bool>::typeName;
}


inline Foam::label Foam::expressions::exprResult::size() const
{
    return size_;
}


template<class Type>
void Foam::expressions::exprResult::setResult
(
    const Field<Type>& val,
    bool wantPointData
)
{
    DebugInFunction << nl;

    target().setResultImpl(val, wantPointData);
}


template<class Type>
void Foam::expressions::exprResult::setResult
(
    Field<Type>&& val,
    bool wantPointData
)
{
    DebugInFunction << nl;

    target().setResultImpl(val, wantPointData);
}


template<class Type>
void Foam::expressions::exprResult::setResultImpl
(
    const Field<Type>& fld,
    bool wantPointData
)
{
    DebugInFunction << nl;

    clear();
    value_.clear();

    isPointData_ = wantPointData;

    size_ = fld.size();
    valType_ = pTraits<Type>::typeName;
    fieldPtr_ = new Field<Type>(fld);

    DebugInFunction << nl;
}


template<class Type>
void Foam::expressions::exprResult::setResultImpl
(
    Field<Type>&& fld,
    bool wantPointData
)
{
    DebugInFunction << nl;

    clear();
    value_.clear();

    isPointData_ = wantPointData;

    size_ = fld.size();
    valType_ = pTraits<Type>::typeName;
    fieldPtr_ = new Field<Type>(std::move(fld));

    DebugInFunction << nl;
}


template<class Type>
void Foam::expressions::exprResult::setResult
(
    Field<Type>* fldPtr,
    bool wantPointData
)
{
    target().setResultImpl(fldPtr, wantPointData);
}


template<class Type>
void Foam::expressions::exprResult::setResultImpl
(
    Field<Type>* fldPtr,
    bool wantPointData
)
{
    clear();
    value_.clear();

    isPointData_ = wantPointData;

    if (fldPtr != nullptr)
    {
        size_ = fldPtr->size();
        valType_ = pTraits<Type>::typeName;
        fieldPtr_ = fldPtr;
    }
}


template<class Type>
void Foam::expressions::exprResult::setResult
(
    const Type& val,
    const label size
)
{
    target().setResultImpl(val, size);
}


template<class Type>
void Foam::expressions::exprResult::setResultImpl
(
    const Type& val,
    const label len
)
{
    DebugInFunction << nl;

    clear();
    value_.set(val);

    isPointData_ = false;

    size_ = len;
    valType_ = pTraits<Type>::typeName;
    fieldPtr_ = new Field<Type>(size_, val);
}


template<class Type>
void Foam::expressions::exprResult::setSingleValue(const Type& val)
{
    target().setSingleValueImpl(val);
}


template<class Type>
bool Foam::expressions::exprResult::writeSingleValueChecked(Ostream& os) const
{
    if (!isType<Type>())
    {
        return false;
    }

    if (this->size() <= 0)
    {
        if (value_.good())
        {
            os << value_.get<Type>();
        }
        else
        {
            // Zero-sized - could write nothing, or a zero value
            os << pTraits<Type>::zero;
        }
    }
    else
    {
        const Field<Type>& fld = *static_cast<const Field<Type>*>(fieldPtr_);

        os << fld.front();
    }

    return true;
}


template<class Type>
bool Foam::expressions::exprResult::writeFieldChecked
(
    const word& keyword,
    Ostream& os
) const
{
    if (!isType<Type>())
    {
        return false;
    }

    if (this->size() <= 0)
    {
        if (value_.good())
        {
            const Type& val = value_.get<Type>();

            if (keyword.empty())
            {
                os << val;
            }
            else
            {
                os.writeEntry(keyword, val);
            }
        }
        else
        {
            // Zero-sized - could write nothing, or a zero value
            if (keyword.empty())
            {
                os << pTraits<Type>::zero;
            }
            else
            {
                Field<Type>().writeEntry(keyword, os);
            }
        }
    }
    else
    {
        const Field<Type>& fld = *static_cast<const Field<Type>*>(fieldPtr_);

        if (keyword.empty())
        {
            os << fld;
        }
        else
        {
            if (value_.good())
            {
                os.writeEntry(keyword, fld.front());
            }
            else
            {
                fld.writeEntry(keyword, os);
            }
        }
    }

    return true;
}


template<class Type>
bool Foam::expressions::exprResult::writeEntryChecked
(
    const word& keyword,
    Ostream& os
) const
{
    if (!isType<Type>())
    {
        return false;
    }

    if (this->size() <= 0)
    {
        if (value_.good() && is_contiguous<Type>::value)
        {
            const Type& val = value_.get<Type>();

            if (keyword.size())
            {
                os.writeKeyword(keyword);
            }
            os  << word("uniform") << token::SPACE << val;
            os.endEntry();
        }
        else
        {
            // Zero-sized - written as nonuniform
            const Field<Type> fld;
            fld.writeEntry(keyword, os);
        }
    }
    else
    {
        const Field<Type>& fld = *static_cast<const Field<Type>*>(fieldPtr_);

        if (value_.good() && is_contiguous<Type>::value)
        {
            if (keyword.size())
            {
                os.writeKeyword(keyword);
            }
            os  << word("uniform") << token::SPACE << fld.front();
            os.endEntry();
        }
        else
        {
            fld.writeEntry(keyword, os);
        }
    }

    return true;
}


template<class Type>
bool Foam::expressions::exprResult::setAverageValueChecked(const bool parRun)
{
    if (!isType<Type>())
    {
        return false;
    }

    const Field<Type>& fld = *static_cast<const Field<Type>*>(fieldPtr_);

    const MinMax<Type> limits = (parRun ? gMinMax(fld) : minMax(fld));

    if (limits.mag() <= SMALL)
    {
        value_.set(limits.centre());  // Average value
    }
    else
    {
        value_.clear();
    }

    return true;
}


template<class Type>
bool Foam::expressions::exprResult::duplicateFieldChecked(const void* ptr)
{
    if (!isType<Type>())
    {
        return false;
    }

    if (fieldPtr_)
    {
        deleteChecked<Type>();
    }

    const Field<Type>& fld = *static_cast<const Field<Type>*>(ptr);

    size_ = fld.size();
    fieldPtr_ = new Field<Type>(fld);

    return true;
}


template<class Type>
void Foam::expressions::exprResult::setSingleValueImpl(const Type& val)
{
    DebugInFunction << nl;

    clear();
    value_.set(val);

    isPointData_ = false;

    size_ = 1;
    valType_ = pTraits<Type>::typeName;
    fieldPtr_ = new Field<Type>(size_, val);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprResult::getResult(bool cacheCopy)
{
    DebugInFunction << nl;

    if (!isType<Type>())
    {
        FatalErrorInFunction
            << "The expected return type " << pTraits<Type>::typeName
            << " is different from the stored result type "
            << valueType() << nl << nl
            << exit(FatalError);
    }

    if (fieldPtr_ == nullptr)
    {
        FatalErrorInFunction
            << "Cannot create tmp from nullptr." << nl
            << "This error message should never appear!!" << nl
            << exit(FatalError);
    }

    Field<Type>* ptr = static_cast<Field<Type>*>(fieldPtr_);

    if (cacheCopy)
    {
        // Leave field intact, return a duplicate field
        // Or return reference instead??
        return tmp<Field<Type>>::New(*ptr);
    }


    tmp<Field<Type>> tresult(ptr);

    fieldPtr_ = nullptr;  // Took ownership of field pointer
    clear();

    return tresult;
}


template<class Type>
inline const Foam::Field<Type>&
Foam::expressions::exprResult::cref() const
{
    DebugInFunction << nl;

    if (!isType<Type>())
    {
        FatalErrorInFunction
            << "The expected return type " << pTraits<Type>::typeName
            << " is different from the stored result type "
            << valueType() << nl << nl
            << exit(FatalError);
    }

    if (fieldPtr_ == nullptr)
    {
        FatalErrorInFunction
            << "Cannot return reference from nullptr." << nl
            << "This error message should never appear!!" << nl
            << exit(FatalError);
    }

    return *static_cast<const Field<Type>*>(fieldPtr_);
}


template<class Type>
inline Foam::Field<Type>&
Foam::expressions::exprResult::ref()
{
    return const_cast<Field<Type>&>(this->cref<Type>());
}


template<class Type>
inline Foam::Field<Type>&
Foam::expressions::exprResult::constCast() const
{
    return const_cast<Field<Type>&>(this->cref<Type>());
}


template<template<class> class BinaryOp, class Type>
inline Type Foam::expressions::exprResult::getReduced
(
    const BinaryOp<Type>& bop,
    const Type& initial
)
{
    if (!isType<Type>())
    {
        FatalErrorInFunction
            << "The expected return type " << pTraits<Type>::typeName
            << " is different from the stored result type "
            << valueType() << nl << nl
            << exit(FatalError);
    }

    Type result = initial;

    const Field<Type>& fld = *static_cast<Field<Type>*>(fieldPtr_);

    for (const Type& val : fld)
    {
        result = bop(result, val);
    }

    return returnReduce(result, bop);
}


// ************************************************************************* //
