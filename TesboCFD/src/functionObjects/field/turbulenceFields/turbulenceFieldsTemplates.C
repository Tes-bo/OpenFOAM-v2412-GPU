/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::turbulenceFields::processField
(
    const word& fieldName,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const word localName(IOobject::scopedName(prefix_, fieldName));

    FieldType* fldPtr = obr_.getObjectPtr<FieldType>(localName);

    if (fldPtr)
    {
        (*fldPtr) == tvalue();
    }
    else
    {
        fldPtr = new FieldType
        (
            IOobject
            (
                localName,
                obr_.time().timeName(),
                obr_,
                IOobject::LAZY_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            tvalue
        );

        obr_.store(fldPtr);
    }
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::nuTilda
(
    const Model& model
) const
{
    const dimensionedScalar omega0(dimless/dimTime, SMALL);

    return volScalarField::New
    (
        "nuTilda.tmp",
        IOobject::NO_REGISTER,
        model.k()/(model.omega() + omega0)
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::L
(
    const Model& model
) const
{
    // (P:Eq. 10.37)
    const scalar Cmu = 0.09;
    const dimensionedScalar eps0(sqr(dimVelocity)/dimTime, SMALL);

    return volScalarField::New
    (
        "L.tmp",
        IOobject::NO_REGISTER,
        pow(Cmu, 0.75)*pow(model.k(), 1.5)/(model.epsilon() + eps0)
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::I
(
    const Model& model
) const
{
    // (P:p. 183)
    // root-mean-square of velocity fluctuations - isotropic turbulence
    const volScalarField uPrime(sqrt((2.0/3.0)*model.k()));
    const dimensionedScalar U0(dimVelocity, SMALL);

    return volScalarField::New
    (
        "I.tmp",
        IOobject::NO_REGISTER,
        uPrime/max(max(uPrime, mag(model.U())), U0)
    );
}


// ************************************************************************* //
