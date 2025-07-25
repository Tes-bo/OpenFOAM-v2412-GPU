/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, f)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    lduInterfaceField(refCast<const lduInterface>(p, dict)),
    faPatchField<Type>(p, iF, dict, requireValue)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coupledFaPatchField<Type>::snGrad() const
{
    return
        (patchNeighbourField() - this->patchInternalField())
       *this->patch().deltaCoeffs();
}


template<class Type>
void Foam::coupledFaPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void Foam::coupledFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    Field<Type>::operator=
    (
        lerp
        (
            this->patchNeighbourField(),
            this->patchInternalField(),
            this->patch().weights()
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -gradientInternalCoeffs();
}


template<class Type>
void Foam::coupledFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
