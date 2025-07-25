/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Henning Scheufler
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "interfaceHeatResistance.H"
#include "constants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"

using namespace Foam::constant;


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::meltingEvaporationModels::
interfaceHeatResistance<Thermo, OtherThermo>
::updateInterface(const volScalarField& T)
{
    const fvMesh& mesh = this->mesh_;

    const volScalarField& alpha = this->pair().from();

    scalarField ap
    (
        volPointInterpolation::New(mesh).interpolate(alpha)
    );

    cutCellIso cutCell(mesh, ap);

    forAll(interfaceArea_, celli)
    {
        label status = cutCell.calcSubCell(celli, isoAlpha_);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] =
                mag(cutCell.faceArea())/mesh.V()[celli];
        }
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp.faceCells(), faceI)
            {
                const label pCelli = pp.faceCells()[faceI];
                bool interface(false);
                if
                (
                     sign(R_.value()) > 0
                 && (T[pCelli] - Tactivate_.value()) > 0
                )
                {
                    interface = true;
                }

                if
                (
                    sign(R_.value()) < 0
                 && (T[pCelli] - Tactivate_.value()) < 0
                )
                {
                    interface = true;
                }

                if
                (
                    interface
                 &&
                    (
                        alpha[pCelli] < 2*isoAlpha_
                     && alpha[pCelli] > 0.5*isoAlpha_
                    )
                )
                {
                    interfaceArea_[pCelli] =
                        mag(pp.faceAreas()[faceI])/mesh.V()[pCelli];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::interfaceHeatResistance<Thermo, OtherThermo>
::interfaceHeatResistance
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    R_("R", dimPower/dimArea/dimTemperature, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    mDotc_
    (
        IOobject
        (
            "mDotc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    mDotcSpread_
    (
        IOobject
        (
            "mDotcSpread",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    htc_
    (
        IOobject
        (
            "htc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimMass/dimArea/dimTemperature/dimTime, Zero)
    ),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5)),
    spread_(dict.getOrDefault<scalar>("spread", 3))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceHeatResistance<Thermo, OtherThermo>
::Kexp(const volScalarField& T)
{

    const fvMesh& mesh = this->mesh_;

    updateInterface(T);

    auto tdeltaT = volScalarField::New
    (
        "tdeltaT",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );
    auto& deltaT = tdeltaT.ref();

    const dimensionedScalar T0(dimTemperature, Zero);

    if (sign(R_.value()) > 0)
    {
        deltaT = max(T - Tactivate_, T0);
    }
    else
    {
        deltaT = max(Tactivate_ - T, T0);
    }

    word fullSpeciesName = this->transferSpecie();
    auto tempOpen = fullSpeciesName.find('.');
    const word speciesName(fullSpeciesName.substr(0, tempOpen));

    tmp<volScalarField> L = mag(this->L(speciesName, T));

    htc_ = R_/L();

    const volScalarField& to = this->pair().to();
    const volScalarField& from = this->pair().from();

    dimensionedScalar D
    (
        "D",
        dimArea,
        spread_/sqr(gAverage(this->mesh_.nonOrthDeltaCoeffs()))
    );

    const dimensionedScalar MdotMin("MdotMin", mDotc_.dimensions(), 1e-3);

    if (max(mDotc_) > MdotMin)
    {
        fvc::spreadSource
        (
            mDotcSpread_,
            mDotc_,
            from,
            to,
            D,
            1e-3
        );
    }

    mDotc_ = interfaceArea_*htc_*deltaT;

    return tmp<volScalarField>::New(mDotc_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceHeatResistance<Thermo, OtherThermo>
::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField coeff(htc_*interfaceArea_);

        if (sign(R_.value()) > 0)
        {
            return(coeff*pos(refValue - Tactivate_));
        }
        else
        {
            return(coeff*pos(Tactivate_ - refValue));
        }
    }

    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceHeatResistance<Thermo, OtherThermo>
::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField coeff(htc_*interfaceArea_*Tactivate_);

        if (sign(R_.value()) > 0)
        {
            return(-coeff*pos(refValue - Tactivate_));
        }
        else
        {
            return(coeff*pos(Tactivate_ - refValue));
        }
    }
    else if (interfaceCompositionModel::P == variable)
    {
        return tmp<volScalarField>::New(mDotcSpread_);
    }

    return nullptr;
}


// ************************************************************************* //
