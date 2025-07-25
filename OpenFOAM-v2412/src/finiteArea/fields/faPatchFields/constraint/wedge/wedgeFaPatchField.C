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

#include "wedgeFaPatch.H"
#include "wedgeFaPatchField.H"
#include "transformField.H"
#include "symmTransform.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    this->evaluate();
}


template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::wedgeFaPatchField<Type>::snGrad() const
{
    const Field<Type> pif (this->patchInternalField());

    return
    (
        transform(refCast<const wedgeFaPatch>(this->patch()).faceT(), pif)
      - pif
    )*(0.5*this->patch().deltaCoeffs());
}


template<class Type>
void Foam::wedgeFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    faPatchField<Type>::operator==
    (
        transform
        (
            refCast<const wedgeFaPatch>(this->patch()).edgeT(),
            this->patchInternalField()
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::wedgeFaPatchField<Type>::snGradTransformDiag() const
{
    const diagTensor diagT =
        0.5*diag(I - refCast<const wedgeFaPatch>(this->patch()).faceT());

    const vector diagV(diagT.xx(), diagT.yy(), diagT.zz());

    return tmp<Field<Type>>::New
    (
        this->size(),
        transformMask<Type>
        (
            pow
            (
                diagV,
                pTraits
                <
                    typename powProduct<vector, pTraits<Type>::rank>::type
                >::zero
            )
        )
    );
}


// ************************************************************************* //
