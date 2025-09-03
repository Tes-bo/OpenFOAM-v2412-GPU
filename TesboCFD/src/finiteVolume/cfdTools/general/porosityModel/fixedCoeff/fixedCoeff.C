/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "fixedCoeff.H"
#include "fvMatrices.H"
#include "pointIndList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porosityModel, fixedCoeff, mesh);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const vectorField& U,
    const scalar rho
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor Cd = rho*(alphaZones[j] + betaZones[j]*mag(U[celli]));
            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);
        }
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorField& AU,
    const vectorField& U,
    const scalar rho
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor alpha = alphaZones[j];
            const tensor beta = betaZones[j];

            AU[celli] += rho*(alpha + beta*mag(U[celli]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const wordRe& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    alphaXYZ_("alpha", dimless/dimTime, coeffs_),
    betaXYZ_("beta", dimless/dimLength, coeffs_),
    alpha_(cellZoneIDs_.size()),
    beta_(cellZoneIDs_.size())
{
    adjustNegativeResistance(alphaXYZ_);
    adjustNegativeResistance(betaXYZ_);

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::calcTransformModelData()
{
    // The alpha coefficient as a tensor
    tensor alphaCoeff(Zero);
    alphaCoeff.xx() = alphaXYZ_.value().x();
    alphaCoeff.yy() = alphaXYZ_.value().y();
    alphaCoeff.zz() = alphaXYZ_.value().z();

    // The beta coefficient as a tensor
    tensor betaCoeff(Zero);
    betaCoeff.xx() = betaXYZ_.value().x();
    betaCoeff.yy() = betaXYZ_.value().y();
    betaCoeff.zz() = betaXYZ_.value().z();

    if (csys().uniform())
    {
        forAll(cellZoneIDs_, zonei)
        {
            alpha_[zonei].resize(1);
            beta_[zonei].resize(1);

            alpha_[zonei] = csys().transform(alphaCoeff);
            beta_[zonei] = csys().transform(betaCoeff);
        }
    }
    else
    {
        forAll(cellZoneIDs_, zonei)
        {
            const pointUIndList cc
            (
                mesh_.cellCentres(),
                mesh_.cellZones()[cellZoneIDs_[zonei]]
            );

            alpha_[zonei] = csys().transform(cc, alphaCoeff);
            beta_[zonei] = csys().transform(cc, betaCoeff);
        }
    }
}


void Foam::porosityModels::fixedCoeff::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), Zero);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();
    const scalar rhoRef = coeffs_.get<scalar>("rhoRef");

    apply(Udiag, Usource, V, U, rhoRef);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.readEntry("rhoRef", rho);
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.readEntry("rhoRef", rho);
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.readEntry("rhoRef", rho);
    }

    apply(AU, U, rho);
}


bool Foam::porosityModels::fixedCoeff::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);

    return true;
}


// ************************************************************************* //
