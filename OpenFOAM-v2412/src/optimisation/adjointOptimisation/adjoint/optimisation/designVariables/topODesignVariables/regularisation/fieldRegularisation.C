/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "fieldRegularisation.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldRegularisation, 1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldRegularisation::fieldRegularisation
(
    fvMesh& mesh,
    const scalarField& alpha,
    const topOZones& zones,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    zones_(zones),
    regularise_(dict.getOrDefault<bool>("regularise", false)),
    project_(dict.getOrDefault<bool>("project", regularise_)),
    radius_(regularisationRadius::New(mesh, dict, false)),
    alpha_(alpha),
    alphaTilda_
    (
        regularise_
      ? new volScalarField
        (
            IOobject
            (
                "alphaTilda",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero),
            fixedValueFvPatchScalarField::typeName
        )
      : nullptr
    ),
    sharpenFunction_
    (
        project_ ?
        topOInterpolationFunction::New(mesh, dict) :
        nullptr
    ),
    regularisationPDE_(regularisationPDE::New(mesh, dict, zones)),
    betaArg_(regularise_ ? alphaTilda_().primitiveField() : alpha_),
    growFromWalls_(dict.getOrDefault<bool>("growFromWalls", false)),
    beta_
    (
        IOobject
        (
            "beta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    )
{
    DebugInfo
        << "Regularise " << Switch(regularise_) << nl
        << "Project " << Switch(project_) << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldRegularisation::updateBeta()
{
    if (regularise_)
    {
        regularise(alpha_, alphaTilda_(), true);
    }

    if (project_)
    {
        sharpenFunction_->interpolate(betaArg_, beta_.primitiveFieldRef());
    }
    else
    {
        beta_.primitiveFieldRef() = betaArg_;
    }

    beta_.correctBoundaryConditions();
}


void Foam::fieldRegularisation::regularise
(
    const scalarField& source,
    scalarField& result,
    const bool isTopoField,
    const regularisationRadius& radius
)
{
    regularisationPDE_->
        regularise(alphaTilda_(), source, result, isTopoField, radius);
}


void Foam::fieldRegularisation::regularise
(
    const scalarField& source,
    scalarField& result,
    const bool isTopoField
)
{
    regularisationPDE_->
        regularise(alphaTilda_(), source, result, isTopoField, radius_());
}


void Foam::fieldRegularisation::postProcessSens(scalarField& sens)
{
    // Add dBeta/dBetaArg
    if (project_)
    {
        sens *= sharpenFunction_->derivative(betaArg_);
    }
    // Add part due to regularisation
    if (regularise_)
    {
        // Solve the adjoint to the regularisation equation
        regularise(sens, sens, false);
    }

    // Add volume
    sens *= mesh_.V();
}


// ************************************************************************* //
