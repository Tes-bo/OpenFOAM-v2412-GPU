/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "thermalShell.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalShell, 0);

addToRunTimeSelectionTable(thermalShellModel, thermalShell, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool thermalShell::init(const dictionary& dict)
{
    if (thickness_ > 0)
    {
        h_ = dimensionedScalar("thickness", dimLength, thickness_);
    }

    this->solution().readEntry("nNonOrthCorr", nNonOrthCorr_);

    return true;
}


tmp<areaScalarField> thermalShell::qr()
{
    auto taqr =
        tmp<areaScalarField>::New
        (
            IOobject
            (
                "tqr",
                regionMesh().time().timeName(),
                regionMesh().thisDb()
            ),
            regionMesh(),
            dimensionedScalar(dimPower/dimArea, Zero)
        );

    if (!qrName_.empty() && qrName_ != "none")
    {
        vsm().mapToSurface<scalar>
        (
            primaryMesh().lookupObject<volScalarField>(qrName_),
            taqr.ref().primitiveFieldRef()
        );
    }

    return taqr;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void thermalShell::solveEnergy()
{
    DebugInFunction << endl;

    const areaScalarField rhoCph(Cp()*rho()*h_);

    faScalarMatrix TEqn
    (
        fam::ddt(rhoCph, T_)
      - fam::laplacian(kappa()*h_, T_)
     ==
        qs_
      + qr()
      + faOptions()(h_, rhoCph, T_)
    );

    TEqn.relax();

    faOptions().constrain(TEqn);

    TEqn.solve();

    faOptions().correct(T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalShell::thermalShell
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    thermalShellModel(modelType, mesh, dict),
    nNonOrthCorr_(1),
    thermo_(dict.subDict("thermo")),
    qs_
    (
        IOobject
        (
            "qs_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPower/dimArea, Zero)
    ),
    h_
    (
        IOobject
        (
            "h_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    thickness_(dict.getOrDefault<scalar>("thickness", 0))
{
    init(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalShell::preEvolveRegion()
{}


void thermalShell::evolveRegion()
{
    nNonOrthCorr_ = solution().get<label>("nNonOrthCorr");

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        solveEnergy();
    }

    Info<< T_.name() << " min/max   = " << gMinMax(T_) << endl;
}


const tmp<areaScalarField> thermalShell::Cp() const
{
    return areaScalarField::New
    (
        "Cps",
        IOobject::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTemperature/dimMass, thermo_.Cp()),
        faPatchFieldBase::zeroGradientType()
    );
}


const tmp<areaScalarField> thermalShell::rho() const
{
    return areaScalarField::New
    (
        "rhos",
        IOobject::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimDensity, thermo_.rho()),
        faPatchFieldBase::zeroGradientType()
    );
}


const tmp<areaScalarField> thermalShell::kappa() const
{
    return areaScalarField::New
    (
        "kappas",
        IOobject::NO_REGISTER,
        regionMesh(),
        dimensionedScalar
        (
            dimPower/dimLength/dimTemperature,
            thermo_.kappa()
        ),
        faPatchFieldBase::zeroGradientType()
    );
}


void thermalShell::info()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
