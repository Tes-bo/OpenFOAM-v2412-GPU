/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "heatExchangerSource.H"
#include "heatExchangerModel.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(heatExchangerSource, 0);
    addToRunTimeSelectionTable(option, heatExchangerSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatExchangerSource::heatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    heatExchangerModelPtr_(heatExchangerModel::New(mesh, name, coeffs_))
{
    read(dict);

    // Set the field name to that of the energy
    // field from which the temperature is obtained

    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();

    heatExchangerModelPtr_->initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::heatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    const scalarField Q(heatExchangerModelPtr_->energyDensity(cells_));

    if (this->V() > VSMALL)
    {
        const word& UName = heatExchangerModelPtr_->U();
        const auto& U = mesh_.lookupObject<volVectorField>(UName);

        const scalarField& V = mesh_.V();
        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            heSource[celli] -= Q[i]*V[celli]*mag(U[celli]);
        }
    }

    heatExchangerModelPtr_->write(log);
}


bool Foam::fv::heatExchangerSource::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict) || !heatExchangerModelPtr_->read(dict))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
