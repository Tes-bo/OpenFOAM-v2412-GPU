/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "cyclicFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{}


template<class Type>
Foam::cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
    if (!isA<cyclicFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    coupledFaPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    cyclicPatch_(refCast<const cyclicFaPatch>(p, dict))
{
    if (!isA<cyclicFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (IOobjectOption::isReadRequired(requireValue))
    {
        this->evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
Foam::cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFaPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& faceCells = cyclicPatch_.faceCells();

    auto tpnf = tmp<Field<Type>>::New(this->size());
    auto& pnf = tpnf.ref();

    const label sizeby2 = this->size()/2;

    if (doTransform())
    {
        for (label facei=0; facei<sizeby2; ++facei)
        {
            pnf[facei] = transform
            (
                forwardT()[0], iField[faceCells[facei + sizeby2]]
            );

            pnf[facei + sizeby2] = transform
            (
                reverseT()[0], iField[faceCells[facei]]
            );
        }
    }
    else
    {
        for (label facei=0; facei<sizeby2; ++facei)
        {
            pnf[facei] = iField[faceCells[facei + sizeby2]];
            pnf[facei + sizeby2] = iField[faceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
void Foam::cyclicFaPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    solveScalarField pnf(this->size());

    const label sizeby2 = this->size()/2;
    const labelUList& faceCells = cyclicPatch_.faceCells();

    for (label facei = 0; facei < sizeby2; ++facei)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFaPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    Field<Type> pnf(this->size());

    const label sizeby2 = this->size()/2;
    const labelUList& faceCells = cyclicPatch_.faceCells();

    for (label facei = 0; facei < sizeby2; ++facei)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


// ************************************************************************* //
