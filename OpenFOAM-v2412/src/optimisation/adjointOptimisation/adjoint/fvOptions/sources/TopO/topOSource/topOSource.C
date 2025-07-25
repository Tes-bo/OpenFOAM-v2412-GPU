/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "topOSource.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "topOVariablesBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(topOSource, 1);
        addToRunTimeSelectionTable
        (
            option,
            topOSource,
            dictionary
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fv::topOSource::getSource()
{
    auto tinterpolant
    (
        tmp<DimensionedField<scalar, volMesh>>::New
        (
            IOobject
            (
                "source",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimless/dimTime,
            scalarField(mesh_.nCells(), Zero)
        )
    );
    DimensionedField<scalar, volMesh>& interpolant = tinterpolant.ref();

    if (mesh_.foundObject<topOVariablesBase>("topOVars"))
    {
        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");
        vars.sourceTerm
            (interpolant, interpolation_(), betaMax_, interpolationFieldName_);

        if (darcyFlow_)
        {
            interpolant.field() += betaMax_*Da_();
        }
    }

    return tinterpolant;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::topOSource::topOSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    interpolation_(topOInterpolationFunction::New(mesh, dict)),
    interpolationFieldName_(word::null),
    betaMax_(0),
    darcyFlow_(false),
    Da_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::topOSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "Adding Brinkman source to " << eqn.psi().name() << endl;

    eqn -= fvm::Sp(getSource(), eqn.psi());
}


void Foam::fv::topOSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "Adding Brinkman source to " << eqn.psi().name() << endl;

    eqn -= fvm::Sp(getSource(), eqn.psi());
}


void Foam::fv::topOSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "Adding Brinkman source to " << eqn.psi().name() << endl;

    eqn -= fvm::Sp(rho*getSource(), eqn.psi());
}


void Foam::fv::topOSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "Adding Brinkman source to " << eqn.psi().name() << endl;

    eqn -= fvm::Sp(rho*getSource(), eqn.psi());
}


void Foam::fv::topOSource::postProcessSens
(
    scalarField& sens,
    const word& fieldName,
    const word& designVariablesName
)
{
    const label fieldi = applyToField(fieldName);
    if
    (
        fieldi != -1
     && mesh_.foundObject<topOVariablesBase>("topOVars")
    )
    {
        DebugInfo
            << "Postprocessing Brinkman sensitivities for field "
            << fieldName << endl;
        const topOVariablesBase& vars =
            mesh_.lookupObject<topOVariablesBase>("topOVars");
        vars.sourceTermSensitivities
        (
            sens,
            interpolation_(),
            betaMax_,
            designVariablesName,
            interpolationFieldName_
         );
    }
}


bool Foam::fv::topOSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        fieldNames_ = coeffs_.get<wordList>("names");
        interpolationFieldName_ = coeffs_.get<word>("interpolationField");
        applied_.setSize(fieldNames_.size(), false);
        if (mesh_.foundObject<topOVariablesBase>("topOVars"))
        {
            const topOVariablesBase& vars =
                mesh_.lookupObject<topOVariablesBase>("topOVars");
            betaMax_ =
                coeffs_.getOrDefault<scalar>("betaMax", vars.getBetaMax());
        }

        darcyFlow_ = coeffs_.getOrDefault<bool>("darcyFlow", false);
        if (darcyFlow_)
        {
            Da_.reset(new scalar(coeffs_.getOrDefault<scalar>("Da", 1.e-5)));
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
