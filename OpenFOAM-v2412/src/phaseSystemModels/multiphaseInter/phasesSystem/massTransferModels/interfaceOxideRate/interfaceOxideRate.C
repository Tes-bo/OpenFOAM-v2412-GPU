/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "interfaceOxideRate.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "timeVaryingMassSorptionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>
::interfaceOxideRate
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_
    (
        dimensionedScalar
        (
            dimDensity/dimTime,
            dict.getCheck<scalar>("C", scalarMinMax::ge(0))
        )
    ),
    Tliquidus_
    (
        dimensionedScalar
        (
            dimTemperature,
            dict.getCheck<scalar>("Tliquidus", scalarMinMax::ge(0))
        )
    ),
    Tsolidus_
    (
        dimensionedScalar
        (
            dimTemperature,
            dict.getCheck<scalar>("Tsolidus", scalarMinMax::ge(0))
        )
    ),
    oxideCrit_
    (
        dimensionedScalar
        (
            dimDensity,
            dict.getCheck<scalar>("oxideCrit", scalarMinMax::ge(0))
        )
    ),
    mDotOxide_
    (
        IOobject
        (
            "mDotOxide",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::Kexp
(
    const volScalarField& T
)
{
    const volScalarField& from = this->pair().from();
    const volScalarField& to = this->pair().to();

    // (CSC:Eq. 2)
    const fvMesh& mesh = this->mesh_;
    scalarField ap
    (
        volPointInterpolation::New(mesh).interpolate(from)
    );

    cutCellIso cutCell(mesh, ap);

    tmp<volScalarField> tSalpha = scalar(0)*from;
    volScalarField& Salpha = tSalpha.ref();

    forAll(Salpha, celli)
    {
        const label status = cutCell.calcSubCell(celli, isoAlpha_);
        if (status == 0) // cell is cut
        {
            Salpha[celli] = scalar(1);
        }
    }

    // (CSC:Eq. 5)
    tmp<volScalarField> tSoxide =
        max((oxideCrit_.value() - to)/oxideCrit_.value(), scalar(0));

    // (CSC:Eq. 4)
    tmp<volScalarField> tST =
        Foam::exp
        (
            scalar(1)
          - scalar(1)/max((T - Tsolidus_)/(Tliquidus_ - Tsolidus_),scalar(1e-6))
        );

    // (CSC:Eq. 6)
    mDotOxide_ = C_*tSalpha*tSoxide*tST;

    const volScalarField::Boundary& alphab = to.boundaryField();

    forAll(alphab, patchi)
    {
        if (isA<timeVaryingMassSorptionFvPatchScalarField>(alphab[patchi]))
        {
            const auto& pp =
                refCast<const timeVaryingMassSorptionFvPatchScalarField>
                (
                    alphab[patchi]
                );
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            tmp<scalarField> tsb = pp.source();

            auto tRhoto = volScalarField::New
            (
                "tRhoto",
                IOobject::NO_REGISTER,
                mesh,
                dimensionedScalar(dimDensity, Zero)
            );
            auto& rhoto = tRhoto.ref();

            rhoto = this->pair().to().rho();

            forAll(fc, faceI)
            {
                const label cellI = fc[faceI];
                const scalar rhoI = rhoto[cellI];
                mDotOxide_[cellI] += rhoI*tsb()[faceI];
            }
        }
    }

    return tmp<volScalarField>::New(mDotOxide_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


// ************************************************************************* //
