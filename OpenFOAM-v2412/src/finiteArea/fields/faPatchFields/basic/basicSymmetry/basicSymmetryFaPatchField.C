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

#include "basicSymmetryFaPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    this->evaluate();
}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFaPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().edgeNormals());

    return
    (
        transform(I - 2.0*sqr(nHat), this->patchInternalField())
      - this->patchInternalField()
    )*(this->patch().deltaCoeffs()/2.0);
}


template<class Type>
void Foam::basicSymmetryFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().edgeNormals());
    Field<Type>::operator=
    (
        (
            this->patchInternalField()
          + transform(I - 2.0*sqr(nHat), this->patchInternalField())
        )/2.0
    );

    transformFaPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFaPatchField<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().edgeNormals()));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// ************************************************************************* //
