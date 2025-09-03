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

#include "PopeIndex.H"
#include "resolutionIndexModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace resolutionIndexModels
{
    defineTypeNameAndDebug(PopeIndex, 0);
    addToRunTimeSelectionTable(resolutionIndexModel, PopeIndex, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::resolutionIndexModels::PopeIndex::kNum() const
{
    const auto& kSgs = getOrReadField<volScalarField>(kName_);
    const auto& Delta = getOrReadField<volScalarField>(deltaName_);

    tmp<volScalarField> th = cbrt(V());

    // (CKJ:Eq. 28)
    return Cn_*sqr(th/Delta)*kSgs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resolutionIndexModels::PopeIndex::PopeIndex
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    resolutionIndexModel(name, mesh, dict),
    includeKnum_(),
    Cn_(),
    UName_(),
    UMeanName_(),
    kName_(),
    deltaName_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::resolutionIndexModels::PopeIndex::read(const dictionary& dict)
{
    if (!resolutionIndexModel::read(dict))
    {
        return false;
    }

    includeKnum_ = dict.getOrDefault<bool>("includeKnum", true);
    if (includeKnum_)
    {
        Cn_ = dict.getOrDefault<scalar>("Cnu", 1.0);
    }
    UName_ = dict.getOrDefault<word>("U", "U");
    UMeanName_ = dict.getOrDefault<word>("UMean", "UMean");
    kName_ = dict.getOrDefault<word>("k", "k");
    deltaName_ = dict.getOrDefault<word>("delta", "delta");

    return true;
}


bool Foam::resolutionIndexModels::PopeIndex::execute()
{
    // Calculate resolved k field
    const auto& U = getOrReadField<volVectorField>(UName_);
    const auto& UMean = getOrReadField<volVectorField>(UMeanName_);
    const volScalarField kRes(0.5*magSqr(U - UMean));

    // Calculate subgrid-scale k field
    const auto& kSgs = getOrReadField<volScalarField>(kName_);

    // Calculate total k field
    tmp<volScalarField> tkTot = kRes + kSgs;
    if (includeKnum_)
    {
        tkTot.ref() += mag(kNum());
    }


    // Calculate index field
    auto& index = getOrReadField<volScalarField>(resultName());

    const dimensionedScalar kMin(kSgs.dimensions(), SMALL);

    // (P:p. 560;CKJ:Eq. 11)
    index = kRes/max(kMin, tkTot);
    index.correctBoundaryConditions();

    return true;
}


bool Foam::resolutionIndexModels::PopeIndex::write()
{
    const auto& index = getOrReadField<volScalarField>(resultName());

    Info<< tab << "writing field:" << index.name() << endl;

    index.write();

    return true;
}


// ************************************************************************* //
