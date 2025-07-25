/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "CelikEtaIndex.H"
#include "resolutionIndexModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace resolutionIndexModels
{
    defineTypeNameAndDebug(CelikEtaIndex, 0);
    addToRunTimeSelectionTable(resolutionIndexModel, CelikEtaIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikEtaIndex::eta() const
{
    const auto& nu = getOrReadField<volScalarField>(nuName_);
    tmp<volScalarField> tepsilon = epsilon();

    const dimensionedScalar epsilonMin(tepsilon.cref().dimensions(), SMALL);

    // (CKJ:Eq. 23)
    return pow025(pow3(nu)/max(epsilonMin, tepsilon));
}


Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikEtaIndex::epsilon() const
{
    const auto& kSgs = getOrReadField<volScalarField>(kName_);
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);
    tmp<volScalarField> tnuEff = nuEff();

    // (Derived based on CKJ:Eq. 25-26, p.031102-5)
    return tnuEff*kSgs/(Ck_*sqr(Delta));
}


Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikEtaIndex::nuEff() const
{
    const auto& nu = getOrReadField<volScalarField>(nuName_);
    const auto& nuSgs = getOrReadField<volScalarField>(nutName_);
    tmp<volScalarField> tnuNum = nuNum();

    // (CKJ:p. 031102-3)
    return tnuNum + nuSgs + nu;
}


Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikEtaIndex::nuNum() const
{
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);

    tmp<volScalarField> tkNum = kNum();

    // (CKJ:Eq. 35)
    return sign(tkNum.cref())*Cnu_*Delta*sqrt(mag(tkNum.cref()));
}


Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikEtaIndex::kNum() const
{
    const auto& kSgs = getOrReadField<volScalarField>(kName_);
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);

    tmp<volScalarField> th = cbrt(V());

    // (CKJ:Eq. 28)
    return Cn_*sqr(th/Delta)*kSgs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resolutionIndexModels::CelikEtaIndex::CelikEtaIndex
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    resolutionIndexModel(name, mesh, dict),
    alphaEta_(),
    m_(),
    Cnu_(),
    Cn_(),
    Ck_(),
    kName_(),
    deltaName_(),
    nuName_(),
    nutName_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::resolutionIndexModels::CelikEtaIndex::read(const dictionary& dict)
{
    if (!resolutionIndexModel::read(dict))
    {
        return false;
    }

    // (Default values from CKJ:p. 031102-3, 031102-5)
    alphaEta_ = dict.getOrDefault<scalar>("alphaEta", 0.05);
    m_ = dict.getOrDefault<scalar>("m", 0.5);
    Cnu_ = dict.getOrDefault<scalar>("Cnu", 0.1);
    Cn_ = dict.getOrDefault<scalar>("Cn", 1.0);
    Ck_ = dict.getOrDefault<scalar>("Ck", 0.0376);
    kName_ = dict.getOrDefault<word>("k", "k");
    deltaName_ = dict.getOrDefault<word>("delta", "delta");
    nuName_ = dict.getOrDefault<word>("nu", "nu");
    nutName_ = dict.getOrDefault<word>("nut", "nut");

    return true;
}


bool Foam::resolutionIndexModels::CelikEtaIndex::execute()
{
    // Calculate Kolmogorov and mesh length scales
    tmp<volScalarField> teta = eta();
    tmp<volScalarField> th = cbrt(V());

    // Calculate index field
    auto& index = getOrReadField<volScalarField>(resultName());

    // (CKJ:Eq. 9)
    index = 1.0/(1.0 + alphaEta_*pow(th/teta, m_));
    index.correctBoundaryConditions();

    return true;
}


bool Foam::resolutionIndexModels::CelikEtaIndex::write()
{
    const auto& index = getOrReadField<volScalarField>(resultName());

    Info<< tab << "writing field:" << index.name() << endl;

    index.write();

    return true;
}


// ************************************************************************* //
