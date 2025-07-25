/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "fixedTemperatureConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fixedTemperatureConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            fixedTemperatureConstraint,
            dictionary
        );
    }
}

const Foam::Enum
<
    Foam::fv::fixedTemperatureConstraint::temperatureMode
>
Foam::fv::fixedTemperatureConstraint::temperatureModeNames_
({
    { temperatureMode::tmUniform, "uniform" },
    { temperatureMode::tmLookup, "lookup" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperatureConstraint::fixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    mode_(temperatureModeNames_.get("mode", coeffs_)),
    Tuniform_(nullptr),
    TName_("T")
{
    switch (mode_)
    {
        case tmUniform:
        {
            Tuniform_.reset
            (
                Function1<scalar>::New("temperature", coeffs_, &mesh_)
            );
            break;
        }
        case tmLookup:
        {
            TName_ = coeffs_.getOrDefault<word>("T", "T");
            break;
        }
        default:
        {
            // Error handling already done by Enum
        }
    }


    // Set the field name to that of the energy
    // field from which the temperature is obtained

    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    switch (mode_)
    {
        case tmUniform:
        {
            const scalar t = mesh_.time().value();
            scalarField Tuni(cells_.size(), Tuniform_->value(t));
            eqn.setValues(cells_, thermo.he(thermo.p(), Tuni, cells_));

            break;
        }
        case tmLookup:
        {
            const auto& T = mesh().lookupObject<volScalarField>(TName_);

            scalarField Tlkp(T, cells_);
            eqn.setValues(cells_, thermo.he(thermo.p(), Tlkp, cells_));

            break;
        }
        default:
        {
            // Error handling already done by Enum
        }
    }
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        if (Tuniform_ && coeffs_.found(Tuniform_->name(), keyType::LITERAL))
        {
            Tuniform_.reset
            (
                Function1<scalar>::New(Tuniform_->name(), dict, &mesh_)
            );
        }

        coeffs_.readIfPresent("T", TName_);

        return true;
    }

    return false;
}


// ************************************************************************* //
