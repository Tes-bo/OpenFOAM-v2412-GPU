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

#include "CelikNuIndex.H"
#include "resolutionIndexModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace resolutionIndexModels
{
    defineTypeNameAndDebug(CelikNuIndex, 0);
    addToRunTimeSelectionTable(resolutionIndexModel, CelikNuIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikNuIndex::nuNum() const
{
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);

    tmp<volScalarField> tkNum = kNum();

    // (CKJ:Eq. 35)
    return sign(tkNum.cref())*Cnu_*Delta*sqrt(tkNum.cref());
}


Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::CelikNuIndex::kNum() const
{
    const auto& kSgs = getOrReadField<volScalarField>(kName_);
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);

    tmp<volScalarField> th = cbrt(V());

    // (CKJ:Eq. 28)
    return Cn_*sqr(th/Delta)*kSgs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resolutionIndexModels::CelikNuIndex::CelikNuIndex
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    resolutionIndexModel(name, mesh, dict),
    alphaNu_(),
    n_(),
    Cnu_(),
    Cn_(),
    kName_(),
    deltaName_(),
    nuName_(),
    nutName_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::resolutionIndexModels::CelikNuIndex::read(const dictionary& dict)
{
    if (!resolutionIndexModel::read(dict))
    {
        return false;
    }

    // (Default values from CKJ:p. 031102-3, 031102-5)
    alphaNu_ = dict.getOrDefault<scalar>("alphaNu", 0.05);
    n_ = dict.getOrDefault<scalar>("n", 0.53);
    Cnu_ = dict.getOrDefault<scalar>("Cnu", 0.1);
    Cn_ = dict.getOrDefault<scalar>("Cn", 1.0);
    kName_ = dict.getOrDefault<word>("k", "k");
    deltaName_ = dict.getOrDefault<word>("delta", "delta");
    nuName_ = dict.getOrDefault<word>("nu", "nu");
    nutName_ = dict.getOrDefault<word>("nut", "nut");

    return true;
}


bool Foam::resolutionIndexModels::CelikNuIndex::execute()
{
    // Calculate effective eddy viscosity field
    const auto& nu = getOrReadField<volScalarField>(nuName_);
    const auto& nuSgs = getOrReadField<volScalarField>(nutName_);
    tmp<volScalarField> tnuNum = nuNum();
    tmp<volScalarField> tnuEff = tnuNum + nuSgs + nu;

    // Calculate index field
    auto& index = getOrReadField<volScalarField>(resultName());

    // (CKJ:Eq. 10)
    index = 1.0/(1.0 + alphaNu_*pow(tnuEff/nu, n_));
    index.correctBoundaryConditions();

    return true;
}


bool Foam::resolutionIndexModels::CelikNuIndex::write()
{
    const auto& index = getOrReadField<volScalarField>(resultName());

    Info<< tab << "writing field:" << index.name() << endl;

    index.write();

    return true;
}


// ************************************************************************* //
