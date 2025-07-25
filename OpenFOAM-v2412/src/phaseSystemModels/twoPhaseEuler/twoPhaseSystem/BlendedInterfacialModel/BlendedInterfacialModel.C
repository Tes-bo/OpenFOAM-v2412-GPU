/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "BlendedInterfacialModel.H"
#include "fixedValueFvsPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class modelType>
template<class GeometricField>
void Foam::BlendedInterfacialModel<modelType>::correctFixedFluxBCs
(
    GeometricField& field
) const
{
    auto& fieldBf = field.boundaryFieldRef();

    forAll(pair_.phase1().phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                pair_.phase1().phi().boundaryField()[patchi]
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::BlendedInterfacialModel
(
    const phasePair::dictTable& modelTable,
    const blendingMethod& blending,
    const phasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1,
    const bool correctFixedFluxBCs
)
:
    pair_(pair),
    pair1In2_(pair1In2),
    pair2In1_(pair2In1),
    blending_(blending),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{
    if (modelTable.found(pair_))
    {
        model_.reset
        (
            modelType::New
            (
                modelTable[pair_],
                pair_
            ).ptr()
        );
    }

    if (modelTable.found(pair1In2_))
    {
        model1In2_.reset
        (
            modelType::New
            (
                modelTable[pair1In2_],
                pair1In2_
            ).ptr()
        );
    }

    if (modelTable.found(pair2In1_))
    {
        model2In1_.reset
        (
            modelType::New
            (
                modelTable[pair2In1_],
                pair2In1_
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<modelType>::K() const
{
    tmp<volScalarField> f1, f2;

    if (model_ || model1In2_)
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_ || model2In1_)
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = volScalarField::New
    (
        IOobject::scopedName(modelType::typeName, "K"),
        IOobject::NO_REGISTER,
        pair_.phase1().mesh(),
        dimensionedScalar(modelType::dimK, Zero)
    );

    if (model_)
    {
        x.ref() += model_->K()*(f1() - f2());
    }

    if (model1In2_)
    {
        x.ref() += model1In2_->K()*(1 - f1);
    }

    if (model2In1_)
    {
        x.ref() += model2In1_->K()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_ || model1In2_ || model2In1_)
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<modelType>::Kf() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_ || model1In2_)
    {
        f1 = fvc::interpolate
        (
            blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    if (model_ || model2In1_)
    {
        f2 = fvc::interpolate
        (
            blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    auto x = surfaceScalarField::New
    (
        IOobject::scopedName(modelType::typeName, "Kf"),
        pair_.phase1().mesh(),
        dimensionedScalar(modelType::dimK, Zero)
    );

    if (model_)
    {
        x.ref() += model_->Kf()*(f1() - f2());
    }

    if (model1In2_)
    {
        x.ref() += model1In2_->Kf()*(1 - f1);
    }

    if (model2In1_)
    {
        x.ref() += model2In1_->Kf()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_ || model1In2_ || model2In1_)
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::BlendedInterfacialModel<modelType>::F() const
{
    tmp<volScalarField> f1, f2;

    if (model_ || model1In2_)
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_ || model2In1_)
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = GeometricField<Type, fvPatchField, volMesh>::New
    (
        IOobject::scopedName(modelType::typeName, "F"),
        IOobject::NO_REGISTER,
        pair_.phase1().mesh(),
        dimensioned<Type>(modelType::dimF, Zero)
    );

    if (model_)
    {
        x.ref() += model_->F()*(f1() - f2());
    }

    if (model1In2_)
    {
        x.ref() += model1In2_->F()*(1 - f1);
    }

    if (model2In1_)
    {
        x.ref() -= model2In1_->F()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_ || model1In2_ || model2In1_)
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<modelType>::Ff() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_ || model1In2_)
    {
        f1 = fvc::interpolate
        (
            blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    if (model_ || model2In1_)
    {
        f2 = fvc::interpolate
        (
            blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    auto x = surfaceScalarField::New
    (
        IOobject::scopedName(modelType::typeName, "Ff"),
        IOobject::NO_REGISTER,
        pair_.phase1().mesh(),
        dimensionedScalar(modelType::dimF*dimArea, Zero)
    );

    x.ref().setOriented();

    if (model_)
    {
        x.ref() += model_->Ff()*(f1() - f2());
    }

    if (model1In2_)
    {
        x.ref() += model1In2_->Ff()*(1 - f1);
    }

    if (model2In1_)
    {
        x.ref() -= model2In1_->Ff()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_ || model1In2_ || model2In1_)
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<modelType>::D() const
{
    tmp<volScalarField> f1, f2;

    if (model_ || model1In2_)
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_ || model2In1_)
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = volScalarField::New
    (
        IOobject::scopedName(modelType::typeName, "D"),
        IOobject::NO_REGISTER,
        pair_.phase1().mesh(),
        dimensionedScalar(modelType::dimD, Zero)
    );

    if (model_)
    {
        x.ref() += model_->D()*(f1() - f2());
    }

    if (model1In2_)
    {
        x.ref() += model1In2_->D()*(1 - f1);
    }

    if (model2In1_)
    {
        x.ref() += model2In1_->D()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_ || model1In2_ || model2In1_)
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
bool Foam::BlendedInterfacialModel<modelType>::hasModel
(
    const class phaseModel& phase
) const
{
    return
    (
        &phase == &(pair_.phase1())
      ? bool(model1In2_)
      : bool(model2In1_)
    );
}


template<class modelType>
const modelType& Foam::BlendedInterfacialModel<modelType>::phaseModel
(
    const class phaseModel& phase
) const
{
    return &phase == &(pair_.phase1()) ? *model1In2_ : *model2In1_;
}


// ************************************************************************* //
