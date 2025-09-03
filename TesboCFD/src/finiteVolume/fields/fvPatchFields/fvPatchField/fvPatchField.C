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

#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fvPatchField<Type>::readValueEntry
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = fvPatchFieldBase::patch();


    const auto* eptr = dict.findEntry("value", keyType::LITERAL);

    if (eptr)
    {
        Field<Type>::assign(*eptr, p.size());
        return true;
    }

    if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Required entry 'value' : missing for patch " << p.name()
            << " in dictionary " << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class Type>
void Foam::fvPatchField<Type>::extrapolateInternal()
{
    fvPatchFieldBase::patch().patchInternalField(internalField_, *this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchFieldBase(p),
    Field<Type>(p.size()),
    internalField_(iF)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& patchType
)
:
    fvPatchFieldBase(p, patchType),
    Field<Type>(p.size()),
    internalField_(iF)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchFieldBase(p),
    Field<Type>(p.size(), value),
    internalField_(iF)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& pfld
)
:
    fvPatchFieldBase(p),
    Field<Type>(pfld),
    internalField_(iF)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    Field<Type>&& pfld
)
:
    fvPatchFieldBase(p),
    Field<Type>(std::move(pfld)),
    internalField_(iF)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    fvPatchFieldBase(p, dict),
    Field<Type>(p.size()),
    internalField_(iF)
{
    readValueEntry(dict, requireValue);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchFieldBase(ptf, p),
    Field<Type>(p.size()),
    internalField_(iF)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        this->extrapolateInternal();
    }
    this->map(ptf, mapper);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf
)
:
    fvPatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(ptf.internalField_)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::check(const fvPatchField<Type>& rhs) const
{
    fvPatchFieldBase::checkPatch(rhs);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatchField<Type>::snGrad() const
{
    return patch().deltaCoeffs()*(*this - patchInternalField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fvPatchField<Type>::patchInternalField() const
{
    return patch().patchInternalField(internalField_);
}


template<class Type>
void Foam::fvPatchField<Type>::patchInternalField(Field<Type>& pfld) const
{
    patch().patchInternalField(internalField_, pfld);
}


template<class Type>
void Foam::fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    Field<Type>& f = *this;

    if (!this->size() && !mapper.distributed())
    {
        f.resize_nocopy(mapper.size());
        if (f.size())
        {
            f = this->patchInternalField();
        }
    }
    else
    {
        // Map all faces provided with mapping data
        Field<Type>::autoMap(mapper);


        // For unmapped faces set to internal field value (zero-gradient)
        if (mapper.hasUnmapped())
        {
            Field<Type> pif(this->patchInternalField());

            if
            (
                mapper.direct()
             && notNull(mapper.directAddressing())
             && mapper.directAddressing().size()
            )
            {
                const labelList& mapAddressing = mapper.directAddressing();

                forAll(mapAddressing, i)
                {
                    if (mapAddressing[i] < 0)
                    {
                        f[i] = pif[i];
                    }
                }
            }
            else if (!mapper.direct() && mapper.addressing().size())
            {
                const labelListList& mapAddressing = mapper.addressing();

                forAll(mapAddressing, i)
                {
                    const labelList& localAddrs = mapAddressing[i];

                    if (!localAddrs.size())
                    {
                        f[i] = pif[i];
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvPatchField<Type>::updateCoeffs()
{
    fvPatchFieldBase::setUpdated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::updateWeightedCoeffs(const scalarField& weights)
{
    // Default behaviour ignores the weights
    if (!updated())
    {
        updateCoeffs();

        fvPatchFieldBase::setUpdated(true);
    }
}


template<class Type>
void Foam::fvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    fvPatchFieldBase::setUpdated(false);
    fvPatchFieldBase::setManipulated(false);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const scalarField& weights
)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label iMatrix,
    const direction cmp
)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (!patchType().empty())
    {
        os.writeEntry("patchType", patchType());
    }
    if (useImplicit())
    {
        os.writeEntry("useImplicit", "true");
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
