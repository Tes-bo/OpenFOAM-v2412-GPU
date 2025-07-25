/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "explicitPorositySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(explicitPorositySource, 0);
    addToRunTimeSelectionTable(option, explicitPorositySource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::explicitPorositySource::explicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    porosityPtr_(nullptr)
{
    read(dict);

    porosityPtr_.reset(porosityModel::New(name_, mesh_, coeffs_, zoneName()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::explicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= alpha*porosityEqn;
}


bool Foam::fv::explicitPorositySource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        if (!coeffs_.readIfPresent("UNames", fieldNames_))
        {
            fieldNames_.resize(1);
            fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
        }

        fv::option::resetApplied();

        return true;
    }

    return false;
}


// ************************************************************************* //
