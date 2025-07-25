/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "variableHeatTransfer.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(variableHeatTransfer, 0);
    addToRunTimeSelectionTable(option, variableHeatTransfer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::variableHeatTransfer::variableHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    UNbrName_(coeffs_.getOrDefault<word>("UNbr", "U")),
    a_(0),
    b_(0),
    c_(0),
    ds_(0),
    Pr_(0),
    AoV_()
{
    if (master_)
    {
        coeffs_.readEntry("a", a_);
        coeffs_.readEntry("b", b_);
        coeffs_.readEntry("c", c_);
        coeffs_.readEntry("ds", ds_);
        coeffs_.readEntry("Pr", Pr_);

        AoV_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "AoV",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::variableHeatTransfer::calculateHtc()
{
    const auto& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    const auto& nbrTurb =
        nbrMesh.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    const auto& nbrThermo =
        nbrMesh.lookupObject<fluidThermo>(basicThermo::dictName);

    const auto& UNbr = nbrMesh.lookupObject<volVectorField>(UNbrName_);

    tmp<volScalarField> ReNbr(mag(UNbr)*ds_*nbrThermo.rho()/nbrTurb.mut());

    tmp<volScalarField> NuNbr(a_*pow(ReNbr, b_)*pow(Pr_, c_));

    scalarField htcNbr(NuNbr*nbrTurb.kappaEff()/ds_);

    tmp<scalarField> htcNbrMapped(interpolate(htcNbr));

    htc_.primitiveFieldRef() = htcNbrMapped*AoV_();
}


bool Foam::fv::variableHeatTransfer::read(const dictionary& dict)
{
    if (interRegionHeatTransferModel::read(dict))
    {
        coeffs_.readIfPresent("UNbr", UNbrName_);

        coeffs_.readIfPresent("a", a_);
        coeffs_.readIfPresent("b", b_);
        coeffs_.readIfPresent("c", c_);
        coeffs_.readIfPresent("ds", ds_);
        coeffs_.readIfPresent("Pr", Pr_);

        return true;
    }

    return false;
}


// ************************************************************************* //
