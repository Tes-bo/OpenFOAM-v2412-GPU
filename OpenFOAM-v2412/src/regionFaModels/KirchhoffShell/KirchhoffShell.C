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

#include "KirchhoffShell.H"
#include "addToRunTimeSelectionTable.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KirchhoffShell, 0);
addToRunTimeSelectionTable(vibrationShellModel, KirchhoffShell, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool KirchhoffShell::init(const dictionary& dict)
{
    this->solution().readEntry("nNonOrthCorr", nNonOrthCorr_);
    return true;
}


void KirchhoffShell::solveDisplacement()
{
    DebugInFunction << endl;

    const Time& time = primaryMesh().time();

    // Create operand fields for solid physics
    const areaScalarField solidMass(rho()*h_);
    const areaScalarField solidD(D()/solidMass);

    // Deep copy old times of shell displacement field
    auto tw0 = tmp<areaScalarField>::New(w_.oldTime());
    auto tw00 = tmp<areaScalarField>::New(w_.oldTime().oldTime());

    // Flag terms to avoid redundant time-derivative calculations
    const bool f0_enabled = (f0_.value() != scalar(0));
    const bool f1_enabled = (f1_.value() != scalar(0));
    const bool f2_enabled = (f2_.value() != scalar(0));

    // Restore various old fields in sub-cycling, if need be
    if (nSubCycles_ > 1)
    {
        w_.oldTime() = w0_;
        w_.oldTime().oldTime() = w00_;

        if (f0_enabled) laplaceW_.oldTime() = laplaceW0_;
        if (f2_enabled) laplace2W_.oldTime() = laplace2W0_;
    }

    for
    (
        subCycleTime wSubCycle
        (
            const_cast<Time&>(time),
            nSubCycles_
        );
       !(++wSubCycle).end();
    )
    {
        laplaceW_ = fac::laplacian(w_);
        laplace2W_ = fac::laplacian(laplaceW_);

        faScalarMatrix wEqn
        (
            fam::d2dt2(w_)
          + solidD*laplace2W_
        ==
            ps_/solidMass
          + faOptions()(solidMass, w_, dimLength/sqr(dimTime))
        );

        // Avoid time-derivative calculations for f terms, if possible
        if (f0_enabled) wEqn -= f0_*sqrt(solidD)*fac::ddt(laplaceW_);
        if (f1_enabled) wEqn += f1_*fam::ddt(w_);
        if (f2_enabled) wEqn += f2_*solidD*fac::ddt(laplace2W_);

        faOptions().constrain(wEqn);

        wEqn.solve();

        // Cache various old fields inside the sub-cycling
        if (wSubCycle.index() >= wSubCycle.nSubCycles())
        {
            w0_ = w_.oldTime();
            w00_ = w_.oldTime().oldTime();

            if (f0_enabled) laplaceW0_ = laplaceW_.oldTime();
            if (f2_enabled) laplace2W0_ = laplace2W_.oldTime();

            // Update shell acceleration
            a_ = fac::d2dt2(w_);
        }
    }

    Info<< w_.name() << " min/max   = " << gMinMax(w_) << endl;

    // Steal the deep-copy of old times to restore the shell displacement
    w_.oldTime() = tw0;
    w_.oldTime().oldTime() = tw00;

    faOptions().correct(w_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KirchhoffShell::KirchhoffShell
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    vibrationShellModel(modelType, mesh, dict),
    ps_
    (
        IOobject
        (
            "ps_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
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
    laplaceW_
    (
        IOobject
        (
            "laplaceW_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(dimLength), Zero)
    ),
    laplace2W_
    (
        IOobject
        (
            "laplace2W_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(pow3(dimLength)), Zero)
    ),
    w0_
    (
        IOobject
        (
            "w0_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    w00_
    (
        IOobject
        (
            "w00_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    laplaceW0_
    (
        IOobject
        (
            "laplaceW0_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(dimLength), Zero)
    ),
    laplace2W0_
    (
        IOobject
        (
            "laplace2W0_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(pow3(dimLength)), Zero)
    ),
    f0_("f0", dimless, dict),
    f1_("f1", inv(dimTime), dict),
    f2_("f2", dimTime, dict),
    nNonOrthCorr_(1),
    nSubCycles_(1)
{
    init(dict);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void KirchhoffShell::preEvolveRegion()
{}


void KirchhoffShell::evolveRegion()
{
    nNonOrthCorr_ = solution().getLabel("nNonOrthCorr");
    nSubCycles_ = solution().getLabel("nSubCycles");

    for (int nonOrth=0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        solveDisplacement();
    }
}


const tmp<areaScalarField> KirchhoffShell::D() const
{
    const dimensionedScalar E("E", dimForce/dimArea , solid().E());
    const dimensionedScalar nu("nu", dimless, solid().nu());

    return tmp<areaScalarField>(E*pow3(h_)/(12*(1 - sqr(nu))));
}


const tmp<areaScalarField> KirchhoffShell::rho() const
{
    return areaScalarField::New
    (
        "rhos",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedScalar("rho", dimDensity, solid().rho()),
        faPatchFieldBase::zeroGradientType()
    );
}

void KirchhoffShell::info()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
