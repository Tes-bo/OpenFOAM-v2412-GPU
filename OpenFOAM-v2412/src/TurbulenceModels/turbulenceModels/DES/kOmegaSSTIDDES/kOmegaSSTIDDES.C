/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "kOmegaSSTIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
const IDDESDelta& kOmegaSSTIDDES<BasicTurbulenceModel>::setDelta() const
{
    if (!isA<IDDESDelta>(this->delta_()))
    {
        FatalErrorInFunction
            << "The delta function must be set to a " << IDDESDelta::typeName
            << " -based model" << exit(FatalError);
    }

    return refCast<const IDDESDelta>(this->delta_());
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::alpha() const
{
    return max(0.25 - this->y_/IDDESDelta_.hmax(), scalar(-5));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::ft
(
    const volScalarField& magGradU
) const
{
    return tanh(pow3(sqr(Ct_)*this->r(this->nut_, magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::fl
(
    const volScalarField& magGradU
) const
{
    return tanh(pow(sqr(Cl_)*this->r(this->nu(), magGradU), 10));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::fdt
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow(Cdt1_*this->r(this->nut_, magGradU), Cdt2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField lRAS(this->lengthScaleRAS());
    const volScalarField lLES(this->lengthScaleLES(CDES));

    const volScalarField alpha(this->alpha());
    const volScalarField expTerm(exp(sqr(alpha)));

    tmp<volScalarField> fB = min(2*pow(expTerm, -9.0), scalar(1));
    const volScalarField fdTilda(max(1 - fdt(magGradU), fB));

    if (fe_)
    {
        tmp<volScalarField> fe1 =
            2*lerp(pow(expTerm, -9.), pow(expTerm, -11.09), pos0(alpha));
        tmp<volScalarField> fe2 = 1 - max(ft(magGradU), fl(magGradU));
        tmp<volScalarField> fe = max(fe1 - 1, scalar(0))*fe2;

        // Original formulation from Shur et al. paper (2008)
        return max
        (
            fdTilda*(1 + fe)*lRAS + (1 - fdTilda)*lLES,
            dimensionedScalar(dimLength, SMALL)
        );
    }

    // Simplified formulation from Gritskevich et al. paper (2011) where fe = 0
    return max
    (
        lerp(lLES, lRAS, fdTilda),
        dimensionedScalar(dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTIDDES<BasicTurbulenceModel>::kOmegaSSTIDDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cdt1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cdt1",
            this->coeffDict_,
            20
        )
    ),
    Cdt2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cdt2",
            this->coeffDict_,
            3
        )
    ),
    Cl_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cl",
            this->coeffDict_,
            5
        )
    ),
    Ct_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct",
            this->coeffDict_,
            1.87
        )
    ),
    fe_
    (
        Switch::getOrAddToDict
        (
            "fe",
            this->coeffDict_,
            true
        )
    ),

    IDDESDelta_(setDelta())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTIDDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        Cdt1_.readIfPresent(this->coeffDict());
        Cdt2_.readIfPresent(this->coeffDict());
        Cl_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::fd() const
{
    const volScalarField alpha(this->alpha());
    const volScalarField expTerm(exp(sqr(alpha)));

    tmp<volScalarField> fB = min(2*pow(expTerm, -9.0), scalar(1));
    return max(1 - fdt(mag(fvc::grad(this->U_))), fB);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
