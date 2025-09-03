/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "proudmanAcousticPower.H"
#include "volFields.H"
#include "basicThermo.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(proudmanAcousticPower, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        proudmanAcousticPower,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::rhoScale
(
    const tmp<volScalarField>& fld
) const
{
    const auto* thermoPtr = getObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        return fld*thermoPtr->rho();
    }

    if (rhoInf_.value() < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << "Incompressible calculation assumed, but no reference density "
            << "set.  Please set the entry 'rhoInf' to an appropriate value"
            << nl
            << exit(FatalError);
    }

    return rhoInf_*fld;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::a() const
{
    const auto* thermoPtr = getObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        const basicThermo& thermo = *thermoPtr;
        return sqrt(thermo.gamma()*thermo.p()/thermo.rho());
    }

    return volScalarField::New
    (
        scopedName("a"),
        IOobject::NO_REGISTER,
        mesh_,
        aRef_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::k() const
{
    if (kName_ != "none")
    {
        return lookupObject<volScalarField>(kName_);
    }

    const auto& turb =
        lookupObject<turbulenceModel>(turbulenceModel::propertiesName);

    return turb.k();
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::epsilon() const
{
    if (epsilonName_ != "none")
    {
        return lookupObject<volScalarField>(epsilonName_);
    }

    if (omegaName_ != "none")
    {
        // Construct epsilon on-the-fly
        const auto& omega = lookupObject<volScalarField>(omegaName_);
        const scalar betaStar = 0.09;
        return betaStar*k()*omega;
    }

    const auto& turb =
        lookupObject<turbulenceModel>(turbulenceModel::propertiesName);

    return turb.epsilon();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::proudmanAcousticPower::proudmanAcousticPower
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    alphaEps_(0.1),
    rhoInf_("0", dimDensity, -1),
    aRef_(dimVelocity, Zero),
    kName_("none"),
    epsilonName_("none"),
    omegaName_("none")
{
    read(dict);

    volScalarField* PAPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName("P_A"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimPower/dimVolume, Zero)
        )
    );

    PAPtr->store();

    volScalarField* LPPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName("L_P"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    LPPtr->store();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::proudmanAcousticPower::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << nl;

        dict.readIfPresent("alphaEps", alphaEps_);
        rhoInf_.readIfPresent("rhoInf", dict);
        aRef_.readIfPresent("aRef", dict);

        if (dict.readIfPresent("k", kName_))
        {
            Info<< "    k field: " << kName_ << endl;
        }
        else
        {
            Info<< "    k field from turbulence model" << endl;
        }

        if (dict.readIfPresent("epsilon", epsilonName_))
        {
            Info<< "    epsilon field: " << epsilonName_ << endl;
        }
        else
        {
            Info<< "    epsilon field from turbulence model (if needed)"
                << endl;
        }

        if (dict.readIfPresent("omega", omegaName_))
        {
            Info<< "    omega field: " << omegaName_ << endl;
        }
        else
        {
            Info<< "    omega field from turbulence model (if needed)" << endl;
        }

        if (epsilonName_ != "none" && omegaName_ != "none")
        {
            FatalIOErrorInFunction(dict)
                << "either epsilon or omega field names can be set but not both"
                << exit(FatalIOError);
        }

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::proudmanAcousticPower::execute()
{
    const volScalarField Mt(sqrt(2*k())/a());

    auto& P_A = mesh_.lookupObjectRef<volScalarField>(scopedName("P_A"));

    P_A = rhoScale(alphaEps_*epsilon()*pow5(Mt));

    auto& L_P = mesh_.lookupObjectRef<volScalarField>(scopedName("L_P"));

    L_P = 10*log10(P_A/dimensionedScalar("PRef", dimPower/dimVolume, 1e-12));

    return true;
}


bool Foam::functionObjects::proudmanAcousticPower::write()
{
    Log << type() << " " << name() << " write:" << nl;

    const auto& P_A = mesh_.lookupObject<volScalarField>(scopedName("P_A"));

    Log << "    writing field " << P_A.name() << nl;

    P_A.write();

    const auto& L_P = mesh_.lookupObject<volScalarField>(scopedName("L_P"));

    Log << "    writing field " << L_P.name() << nl;

    L_P.write();

    Log << endl;

    return true;
}


// ************************************************************************* //
