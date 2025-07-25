/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020,2023 OpenCFD Ltd
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

#include "buoyancyTurbSource.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), SMALL);

    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tnut(turbPtr->nut());
    const volScalarField& nut = tnut();

    const dictionary& turbDict = turbPtr->coeffDict();
    const scalar Prt
    (
        turbDict.getCheckOrDefault<scalar>("Prt", 0.85, scalarMinMax::ge(SMALL))
    );

    // (DTR:Eq. 21, buoyancy correction term)
    const tmp<volScalarField> GbByK((nut/Prt)*(fvc::grad(rho) & g_)/max(k, k0));

    eqn -= fvm::Sp(GbByK, k);
}


// ************************************************************************* //
