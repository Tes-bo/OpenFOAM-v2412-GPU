/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "PoissonPatchDistMethod.H"
#include "fvcGrad.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(Poisson, 0);
    addToRunTimeSelectionTable(patchDistMethod, Poisson, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::Poisson::Poisson
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


Foam::patchDistMethods::Poisson::Poisson
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::Poisson::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::Poisson::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    if (!tyPsi_)
    {
        tyPsi_ = volScalarField::New
        (
            "yPsi",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(sqr(dimLength), Zero),
            y.boundaryFieldRef().types()
        );
    }
    auto& yPsi = tyPsi_.ref();

    solve(fvm::laplacian(yPsi) == dimensionedScalar("1", dimless, -1.0));

    volVectorField gradyPsi(fvc::grad(yPsi));
    volScalarField magGradyPsi(mag(gradyPsi));

    // Need to stabilise the y for overset meshes since the holed cells
    // keep the initial value (0.0) so the gradient of that will be
    // zero as well. Turbulence models do not like zero wall distance.
    y = max
    (
        sqrt(magSqr(gradyPsi) + 2*yPsi) - magGradyPsi,
        dimensionedScalar("smallY", dimLength, SMALL)
    );

    // For overset: enforce smooth y field (yPsi is smooth, magGradyPsi is
    // not)
    mesh_.interpolate(y);

    // Need to stabilise the y for overset meshes since the holed cells
    // keep the initial value (0.0) so the gradient of that will be
    // zero as well. Turbulence models do not like zero wall distance.
    y.clamp_min(SMALL);

    // For overset: enforce smooth y field (yPsi is smooth, magGradyPsi is
    // not)
    mesh_.interpolate(y);

    // Cache yPsi if the mesh is moving otherwise delete
    if (!mesh_.changing())
    {
        tyPsi_.clear();
    }

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n =
           -gradyPsi
           /max
            (
                magGradyPsi,
                dimensionedScalar("smallMagGradyPsi", dimLength, SMALL)

            );

        // For overset: enforce smooth field
        mesh_.interpolate(n);
    }

    return true;
}


// ************************************************************************* //
