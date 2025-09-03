/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "hydrostaticPressure.H"
#include "basicThermo.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "constrainPressure.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(hydrostaticPressure, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        hydrostaticPressure,
        dictionary
    );
}
}


Foam::dimensionedScalar
Foam::functionObjects::hydrostaticPressure::pRef() const
{
    if (pRefName_ == "none")
    {
        return dimensionedScalar(dimPressure, Zero);
    }
    else if (pRefName_ == "pInf")
    {
        return dimensionedScalar("pRef", dimPressure, pRefValue_);
    }
    else
    {
        return mesh_.lookupObject<uniformDimensionedScalarField>(pRefName_);
    }
}


void Foam::functionObjects::hydrostaticPressure::calculateAndWrite()
{
    const auto& pRef = this->pRef();
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& gh = mesh_.lookupObject<volScalarField>(ghName_);
    const auto& ghf = mesh_.lookupObject<surfaceScalarField>(ghfName_);
    auto& rho = mesh_.lookupObjectRef<volScalarField>(rhoName_);
    auto& thermo = mesh_.lookupObjectRef<basicThermo>(basicThermo::dictName);
    auto& p_rgh = mesh_.lookupObjectRef<volScalarField>(p_rghName_);
    auto& ph_rgh = mesh_.lookupObjectRef<volScalarField>(ph_rghName_);

    auto& p = thermo.p();

    Info<< "Performing hydrostatic pressure initialisation";
    if (!mesh_.regionName().empty())
    {
        Info<< " region=" << mesh_.name();
    }


    if (thermo.incompressible())
    {
        Info<< ": incompressible" << endl;

        // Constant density and temperature

        thermo.correct();
        rho = thermo.rho();
        p = ph_rgh + rho*gh + pRef;
        p_rgh = ph_rgh;
    }
    else
    {
        Info<< ": compressible" << endl;

        p = ph_rgh + rho*gh + pRef;
        thermo.correct();
        rho = thermo.rho();

        for (label i=0; i<nCorrectors_; ++i)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));

            surfaceScalarField phig
            (
                "phig",
               -rhof*ghf*fvc::snGrad(rho)*mesh_.magSf()
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
            );

            ph_rghEqn.relax();

            ph_rghEqn.solve();

            p = ph_rgh + rho*gh + pRef;
            thermo.correct();
            rho = thermo.rho();

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;
        }

        p_rgh = ph_rgh;

        Log << "    writing field " << ph_rgh.name() << nl;
        ph_rgh.write();
    }

    Log << "    writing field " << rho.name() << nl;
    rho.write();

    Log << "    writing field " << p_rgh.name() << nl;
    p_rgh.write();

    Log << "    writing field " << p.name() << nl;
    p.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::hydrostaticPressure::hydrostaticPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    p_rghName_("p_rgh"),
    ph_rghName_("ph_rgh"),
    pRefName_("pRef"),
    pRefValue_(0),
    rhoName_("rho"),
    UName_("U"),
    ghName_("gh"),
    ghfName_("ghf"),
    nCorrectors_(5)
{
    if (read(dict))
    {
        // Read and store the initial ph_rgh field
        volScalarField* ptr =
            new volScalarField
            (
                IOobject
                (
                    ph_rghName_,
                    runTime.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,  // To enable restart
                    IOobject::REGISTER
                ),
                mesh_
            );

        regIOobject::store(ptr);

        bool reInitialise = dict.getOrDefault("reInitialise", false);

        if (runTime.timeIndex() == 0 || reInitialise)
        {
            calculateAndWrite();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::hydrostaticPressure::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readIfPresent("p_rgh", p_rghName_);
        dict.readIfPresent("ph_rgh", ph_rghName_);
        dict.readIfPresent("pRef", pRefName_);
        dict.readIfPresent("rho", rhoName_);
        dict.readIfPresent("U", UName_);
        dict.readIfPresent("gh", ghName_);
        dict.readIfPresent("ghf", ghfName_);
        dict.readIfPresent("nCorrectors", nCorrectors_);

        pRefValue_ = 0;
        if (pRefName_ == "pInf")
        {
            pRefValue_ = dict.get<scalar>("pRefValue");
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::hydrostaticPressure::execute()
{
    return true;
}


bool Foam::functionObjects::hydrostaticPressure::write()
{
    return true;
}


// ************************************************************************* //
