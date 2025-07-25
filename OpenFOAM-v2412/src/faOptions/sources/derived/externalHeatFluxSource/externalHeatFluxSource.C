/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "externalHeatFluxSource.H"
#include "fam.H"
#include "faScalarMatrix.H"
#include "physicoChemicalConstants.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(externalHeatFluxSource, 0);
    addToRunTimeSelectionTable(option, externalHeatFluxSource, dictionary);
}
}


const Foam::Enum
<
    Foam::fa::externalHeatFluxSource::operationMode
>
Foam::fa::externalHeatFluxSource::operationModeNames
({
    { operationMode::fixedPower, "power" },
    { operationMode::fixedHeatFlux, "flux" },
    { operationMode::fixedHeatTransferCoeff, "coefficient" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::externalHeatFluxSource::externalHeatFluxSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& m
)
:
    fa::faceSetOption(sourceName, modelType, dict, m),
    mode_(operationModeNames.get("mode", dict)),
    TName_(dict.getOrDefault<word>("T", "T")),
    Q_(nullptr),
    q_(nullptr),
    h_(nullptr),
    Ta_(nullptr),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 0))
{
    fieldNames_.resize(1, TName_);

    fa::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::externalHeatFluxSource::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo
            << name() << ": applying source to "
            << eqn.psi().name() << endl;

        scalar qflux = 0;

        const scalar timeVal = mesh_.time().timeOutputValue();

        switch (mode_)
        {
            case fixedPower:
            {
                // From [W] to [W/m2]
                qflux = Q_->value(timeVal)/(faceSetOption::A() + VSMALL);
                break;
            }

            case fixedHeatFlux:
            {
                qflux = q_->value(timeVal);
                break;
            }

            default:
            {
                break;
            }
        }

        switch (mode_)
        {
            case fixedPower:
            case fixedHeatFlux:
            {
                auto tQ = DimensionedField<scalar, areaMesh>::New
                (
                    "Q",
                    regionMesh(),
                    dimensionedScalar(dimPower/sqr(dimLength), Zero)
                );
                auto& Q = tQ.ref();

                if (faceSetOption::useSubMesh())
                {
                    UIndirectList<scalar>(Q.field(), faceSetOption::faces())
                        = qflux;
                }
                else
                {
                    Q.field() = qflux;
                }

                eqn += Q;

                break;
            }

            case fixedHeatTransferCoeff:
            {
                const dimensionedScalar Ta
                (
                    "Ta",
                    dimTemperature,
                    Ta_->value(timeVal)
                );

                auto thp = DimensionedField<scalar, areaMesh>::New
                (
                    "h",
                    regionMesh(),
                    dimensionedScalar
                    (
                        "h",
                        dimPower/sqr(dimLength)/dimTemperature,
                        h_->value(timeVal)
                    )
                );
                auto& hp = thp.ref();

                DimensionedField<scalar, areaMesh> hpTa(hp*Ta);

                if (emissivity_ > 0)
                {
                    hp -= emissivity_*sigma.value()*pow3(eqn.psi());
                }

                // Zero htc for non-mapped faces
                faceSetOption::subsetFilter(hp.field());
                faceSetOption::subsetFilter(hpTa.field());

                eqn -= fam::SuSp(hp, eqn.psi()) - hpTa;
                break;
            }
        }
    }
}


bool Foam::fa::externalHeatFluxSource::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        dict.readIfPresent("T", TName_);
        dict.readIfPresent("emissivity", emissivity_);

        mode_ = operationModeNames.get("mode", dict);

        switch (mode_)
        {
            case fixedPower:
            {
                Q_ = Function1<scalar>::New("Q", dict, &mesh_);
                break;
            }
            case fixedHeatFlux:
            {
                Q_ = Function1<scalar>::New("q", dict, &mesh_);
                break;
            }
            case fixedHeatTransferCoeff:
            {
                h_ = Function1<scalar>::New("h", dict, &mesh_);
                Ta_ = Function1<scalar>::New("Ta", dict, &mesh_);
                break;
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
