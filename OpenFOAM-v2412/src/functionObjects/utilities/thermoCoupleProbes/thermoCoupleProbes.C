/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "thermoCoupleProbes.H"
#include "mathematicalConstants.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermoCoupleProbes, 0);
    addToRunTimeSelectionTable(functionObject, thermoCoupleProbes, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermoCoupleProbes::thermoCoupleProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    probes(name, runTime, dict, loadFromFiles, false),
    ODESystem(),
    UName_(dict.getOrDefault<word>("U", "U")),
    radiationFieldName_(dict.get<word>("radiationField")),
    thermo_(mesh_.lookupObject<fluidThermo>(basicThermo::dictName)),
    odeSolver_(nullptr),
    Ttc_()
{
    if (readFields)
    {
        read(dict);
    }

    // Check if the property exists (resume old calculation) or is new
    dictionary probeDict;
    if (getDict(typeName, probeDict))
    {
        probeDict.readEntry("Tc", Ttc_);
    }
    else
    {
        Ttc_ = probes::sample(thermo_.T());
    }

    // Note: can only create the solver once all samples have been found
    // - the number of samples is used to set the size of the ODE system
    odeSolver_ = ODESolver::New(*this, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::label Foam::functionObjects::thermoCoupleProbes::nEqns() const
{
    return this->size();
}


void Foam::functionObjects::thermoCoupleProbes::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    scalarField G(y.size(), Zero);
    scalarField Tc(y.size(), Zero);
    scalarField Uc(y.size(), Zero);
    scalarField rhoc(y.size(), Zero);
    scalarField muc(y.size(), Zero);
    scalarField Cpc(y.size(), Zero);
    scalarField kappac(y.size(), Zero);

    if (radiationFieldName_ != "none")
    {
        G = sample(mesh_.lookupObject<volScalarField>(radiationFieldName_));
    }

    Tc = probes::sample(thermo_.T());

    Uc = mag(this->sample(mesh_.lookupObject<volVectorField>(UName_)));

    rhoc = this->sample(thermo_.rho()());
    kappac = this->sample(thermo_.kappa()());
    muc = this->sample(thermo_.mu()());
    Cpc = this->sample(thermo_.Cp()());

    scalarField Re(rhoc*Uc*d_/muc);
    scalarField Pr(Cpc*muc/kappac);
    Pr = max(ROOTVSMALL, Pr);
    scalarField Nu(2.0 + (0.4*sqrt(Re) + 0.06*pow(Re, 2.0/3.0))*pow(Pr, 0.4));
    scalarField htc(Nu*kappac/d_);

    const scalar sigma = physicoChemical::sigma.value();

    scalar area = 4*constant::mathematical::pi*sqr(0.5*d_);
    scalar volume = (4.0/3.0)*constant::mathematical::pi*pow3(0.5*d_);

    dydx =
        (epsilon_*(G/4 - sigma*pow4(y))*area + htc*(Tc - y)*area)
      / (rho_*Cp_*volume);
}


void Foam::functionObjects::thermoCoupleProbes::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    derivatives(x, y, dfdt);

    const label n = nEqns();

    for (label i=0; i<n; i++)
    {
        for (label j=0; j<n; j++)
        {
            dfdy(i, j) = 0.0;
        }
    }
}


bool Foam::functionObjects::thermoCoupleProbes::write()
{
    if (!pointField::empty())
    {
        (void) prepare(ACTION_WRITE);

        const auto& Tfield = thermo_.T();
        writeValues(Tfield.name(), Ttc_, Tfield.time().timeOutputValue());

        dictionary probeDict;
        probeDict.add("Tc", Ttc_);
        setProperty(typeName, probeDict);
        return true;
    }

    return false;
}


bool Foam::functionObjects::thermoCoupleProbes::execute()
{
    if (!pointField::empty())
    {
        scalar dt = mesh_.time().deltaTValue();
        scalar t = mesh_.time().value();
        odeSolver_->solve(t, t+dt, Ttc_, dt);
        return true;
    }

    return false;
}


bool Foam::functionObjects::thermoCoupleProbes::read(const dictionary& dict)
{
    if (probes::read(dict))
    {
        dict.readEntry("rho", rho_);
        dict.readEntry("Cp", Cp_);
        dict.readEntry("d", d_);
        dict.readEntry("epsilon", epsilon_);
        return true;
    }

    return false;
}


// ************************************************************************* //
