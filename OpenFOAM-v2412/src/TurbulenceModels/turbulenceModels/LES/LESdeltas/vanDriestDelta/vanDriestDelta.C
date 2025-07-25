/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "vanDriestDelta.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDistAddressing.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(vanDriestDelta, 0);
    addToRunTimeSelectionTable(LESdelta, vanDriestDelta, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::vanDriestDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U = turbulenceModel_.U();
    const tmp<volScalarField> tnu = turbulenceModel_.nu();
    const volScalarField& nu = tnu();
    tmp<volScalarField> nuSgs = turbulenceModel_.nut();

    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh.time().constant(),
            mesh.thisDb(),
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("ystar", dimLength, GREAT)
    );

    const fvPatchList& patches = mesh.boundary();
    volScalarField::Boundary& ystarBf = ystar.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U.boundaryField()[patchi];
            const scalarField& nuw = nu.boundaryField()[patchi];
            const scalarField& nuSgsw = nuSgs().boundaryField()[patchi];

            ystarBf[patchi] =
                nuw/sqrt((nuw + nuSgsw)*mag(Uw.snGrad()) + VSMALL);
        }
    }

    // Construct wall transporter
    const auto& wDist = wallDistAddressing::New(mesh);

    // Get distance to nearest wall
    const auto& y = wDist.y();

    // Get ystar from nearest wall
    wDist.map(ystar, mapDistribute::transform());

    // Calculate y/ystar (stored in ystar!) and do the clipping
    constexpr scalar yPlusCutOff = 500;
    // Allow for some precision loss from transformation/interpolation of GREAT
    // (= unvisited value)(though ystar is scalar so should not be transformed)
    constexpr scalar fuzzyGREAT = 0.5*GREAT;

    ystar.dimensions().reset(y.dimensions()/ystar.dimensions());
    forAll(y, celli)
    {
        const scalar yPlus = y[celli]/ystar[celli];
        if (y[celli] > fuzzyGREAT || (yPlus > yPlusCutOff))
        {
            // unvisited : y is GREAT, ystar is 1.0
            ystar[celli] = GREAT;
        }
        else
        {
            ystar[celli] = yPlus;
        }
    }

    forAll(y.boundaryField(), patchi)
    {
        const auto& yp = y.boundaryField()[patchi];
        auto& ystarp = ystar.boundaryFieldRef()[patchi];

        forAll(yp, i)
        {
            const scalar yPlus = yp[i]/ystarp[i];
            if (yp[i] > fuzzyGREAT || (yPlus > yPlusCutOff))
            {
                ystarp[i] = GREAT;
            }
            else
            {
                ystarp[i] = yPlus;
            }
        }
    }
    ystar.correctBoundaryConditions();

    // Note: y/ystar stored in ystar
    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        //(kappa_/Cdelta_)*((scalar(1) + SMALL) - exp(-y/ystar/Aplus_))*y
        (kappa_/Cdelta_)*((scalar(1) + SMALL) - exp(-ystar/Aplus_))*y
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::vanDriestDelta::vanDriestDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    geometricDelta_
    (
        LESdelta::New
        (
            IOobject::groupName("geometricDelta", turbulence.U().group()),
            turbulence,
            // Note: cannot use optionalSubDict - if no *Coeffs dict the
            // code will get stuck in a loop attempting to read the delta entry
            // - consider looking up "geometricDelta" instead of "delta"?
            dict.subDict(type() + "Coeffs")
        )
    ),
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    Aplus_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Aplus",
            26.0
        )
    ),
    Cdelta_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Cdelta",
            0.158
        )
    )
{
    calcInterval_ = 1;
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");
    if (!coeffsDict.readIfPresent("calcInterval", calcInterval_))
    {
        coeffsDict.readIfPresent("updateInterval", calcInterval_);
    }

    delta_ = geometricDelta_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::vanDriestDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");

    geometricDelta_().read(coeffsDict);
    dict.readIfPresent<scalar>("kappa", kappa_);
    coeffsDict.readIfPresent<scalar>("Aplus", Aplus_);
    coeffsDict.readIfPresent<scalar>("Cdelta", Cdelta_);
    calcInterval_ = 1;
    if (!coeffsDict.readIfPresent<label>("calcInterval", calcInterval_))
    {
        coeffsDict.readIfPresent("updateInterval", calcInterval_);
    }

    calcDelta();
}


void Foam::LESModels::vanDriestDelta::correct()
{
    if
    (
        (calcInterval_ > 0)
     && (turbulenceModel_.mesh().time().timeIndex() % calcInterval_) == 0
    )
    {
        geometricDelta_().correct();
        calcDelta();
    }
}


// ************************************************************************* //
