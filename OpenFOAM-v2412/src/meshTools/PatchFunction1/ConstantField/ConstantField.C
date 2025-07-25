/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "ConstantField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& entryName,
    const Type& uniformValue,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    isUniform_(true),
    uniformValue_(uniformValue),
    value_(this->size(), uniformValue_)
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& entryName,
    const bool isUniform,
    const Type& uniformValue,
    const Field<Type>& fieldValues,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    isUniform_(isUniform),
    uniformValue_(uniformValue),
    value_(fieldValues)
{
    if (fieldValues.size() != this->size())
    {
        FatalIOErrorInFunction(dict)
            << "Supplied field size " << fieldValues.size()
            << " is not equal to the number of "
            << (faceValues ? "faces" : "points") << ' '
            << this->size() << " of patch " << pp.name() << nl
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::Field<Type>
Foam::PatchFunction1Types::ConstantField<Type>::getValue
(
    const word& entryName,
    const entry* eptr,
    const dictionary& dict,
    const label len,
    bool& isUniform,
    Type& uniformValue
)
{
    isUniform = true;
    uniformValue = Zero;

    Field<Type> fld;

    if (!eptr)
    {
        if (entryName == dict.dictName())
        {
            // The containing dictionary *IS* the Function1 entry
            // Eg,
            // uniformValue { type constant; value 1.2; }

            dict.readEntry("value", uniformValue);

            fld.resize(len);
            fld = uniformValue;
            return fld;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Null entry" << nl
                << exit(FatalIOError);
        }
    }
    else if (!eptr->isStream())
    {
        // Dictionary format. Eg,
        // key { type constant; value 1.2; }

        dict.readEntry("value", uniformValue);

        fld.resize(len);
        fld = uniformValue;
        return fld;
    }

    ITstream& is = eptr->stream();

    if (is.peek().isWord())
    {
        const word contentType(is);

        if (contentType == "constant" || contentType == "uniform")
        {
            is >> uniformValue;
            fld.resize(len);
            fld = uniformValue;
        }
        else if (contentType == "nonuniform")
        {
            if (len)
            {
                isUniform = false;
            }

            is >> static_cast<List<Type>&>(fld);
            const label lenRead = fld.size();
            if (len != lenRead)
            {
                if
                (
                    len < lenRead
                 && FieldBase::allowConstructFromLargerSize
                )
                {
                    #ifdef FULLDEBUG
                    IOWarningInFunction(dict)
                        << "Sizes do not match. Truncating " << lenRead
                        << " entries to " << len << endl;
                    #endif

                    // Truncate the data
                    fld.resize(len);
                }
                else
                {
                    FatalIOErrorInFunction(dict)
                        << "size " << lenRead
                        << " is not equal to the expected length " << len
                        << exit(FatalIOError);
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Expected keyword 'constant', 'uniform', or 'nonuniform'"
                << ", found " << contentType
                << exit(FatalIOError);
        }
    }
    else
    {
        // Uniform (constant) field
        is >> uniformValue;
        fld.resize(len);
        fld = uniformValue;
    }

    return fld;
}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    isUniform_(true),  // overwritten by getValue()
    uniformValue_(Zero),  // overwritten by getValue()
    value_
    (
        getValue
        (
            entryName,
            dict.findEntry(entryName, keyType::LITERAL),
            dict,
            this->size(),
            isUniform_,
            uniformValue_
        )
    )
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const entry* eptr,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    isUniform_(true),  // overwritten by getValue()
    uniformValue_(Zero),  // overwritten by getValue()
    value_
    (
        getValue
        (
            entryName,
            eptr,
            dict,
            this->size(),
            isUniform_,
            uniformValue_
        )
    )
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const ConstantField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    isUniform_(rhs.isUniform_),
    uniformValue_(rhs.uniformValue_),
    value_(rhs.value_)
{
    // If sizes are different...
    value_.resize(this->size(), Zero);

    if (isUniform_)
    {
        value_ = uniformValue_;
    }
}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const ConstantField<Type>& rhs
)
:
    ConstantField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    value_.autoMap(mapper);

    // If originating from single value override just to make sure
    if (isUniform_)
    {
        value_ = uniformValue_;
    }
}


template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    const auto& cst = refCast<const ConstantField<Type>>(pf1);
    value_.rmap(cst.value_, addr);
}


template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    if (isUniform_)
    {
        os.writeKeyword(this->name())
            << word("constant") << token::SPACE << uniformValue_;
        os.endEntry();
    }
    else
    {
        value_.writeEntry(this->name(), os);
    }
}


// ************************************************************************* //
