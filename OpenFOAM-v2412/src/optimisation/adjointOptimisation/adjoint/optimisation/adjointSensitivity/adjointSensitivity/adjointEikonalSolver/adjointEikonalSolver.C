/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "adjointEikonalSolver.H"
#include "adjointSolver.H"
#include "fvc.H"
#include "fvm.H"
#include "surfaceInterpolation.H"
#include "volFieldsFwd.H"
#include "wallFvPatch.H"
#include "patchDistMethod.H"
#include "fvOptions.H"
#include "zeroGradientFvPatchField.H"
#include "sensitivityTopO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointEikonalSolver, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

wordList adjointEikonalSolver::patchTypes() const
{
    wordList daTypes
    (
        mesh_.boundary().size(),
        fixedValueFvPatchScalarField::typeName
    );

    for (const label patchi : wallPatchIDs_)
    {
        daTypes[patchi] = fvPatchFieldBase::zeroGradientType();
    }

    return daTypes;
}


void adjointEikonalSolver::read()
{
    nEikonalIters_ = dict_.getOrDefault<label>("iters", 1000);
    tolerance_ = dict_.getOrDefault<scalar>("tolerance", 1e-6);
    const scalar defaultEps =
        mesh_.schemesDict().subDict("wallDist").
            subOrEmptyDict("advectionDiffusionCoeffs").
                getOrDefault<scalar>("epsilon", 0.1);
    epsilon_ = dict_.getOrDefault<scalar>("epsilon", defaultEps);
}


tmp<surfaceScalarField> adjointEikonalSolver::computeYPhi()
{
    // Primal distance field
    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();

    volVectorField ny
    (
        IOobject
        (
            "ny",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        patchDistMethod::patchTypes<vector>(mesh_, wallPatchIDs_)
    );

    const fvPatchList& patches = mesh_.boundary();
    volVectorField::Boundary& nybf = ny.boundaryFieldRef();

    for (const label patchi : wallPatchIDs_)
    {
        nybf[patchi] == -patches[patchi].nf();
    }

    ny = fvc::grad(d);

    surfaceVectorField nf(fvc::interpolate(ny));

    return tmp<surfaceScalarField>::New("yPhi", mesh_.Sf() & nf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointEikonalSolver::adjointEikonalSolver
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver,
    const labelHashSet& sensitivityPatchIDs
)
:
    mesh_(mesh),
    dict_(dict.subOrEmptyDict("adjointEikonalSolver")),
    adjointSolver_(adjointSolver),
    sensitivityPatchIDs_(sensitivityPatchIDs),
    nEikonalIters_(-1),
    tolerance_(-1),
    epsilon_(Zero),
    wallPatchIDs_(mesh_.boundaryMesh().findPatchIDs<wallPolyPatch>()),
    da_
    (
        IOobject
        (
            word
            (
                adjointSolver.useSolverNameForFields() ?
                "da" + adjointSolver.solverName() :
                "da"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
//            adjointVars.writeFields() ?
//              IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(adjointSolver.daDimensions() , Zero),
        patchTypes()
    ),
    source_
    (
        IOobject
        (
            "sourceEikonal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(adjointSolver.daDimensions()/dimLength, Zero)
    ),
    distanceSensPtr_(createZeroBoundaryPtr<vector>(mesh_))
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool adjointEikonalSolver::readDict(const dictionary& dict)
{
    dict_ = dict.subOrEmptyDict("adjointEikonalSolver");
    read();

    return true;
}


void adjointEikonalSolver::accumulateIntegrand(const scalar dt)
{
    // Accumulate integrand from the current time step
    source_ += adjointSolver_.adjointEikonalSource()*dt;
}


void adjointEikonalSolver::solve()
{
    read();

    // Primal distance field
    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();

    // Convecting flux
    tmp<surfaceScalarField> tyPhi = computeYPhi();
    const surfaceScalarField& yPhi = tyPhi();

    volScalarField scaleDims
    (
        IOobject
        (
            "scaleDims",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar("scaleDims", dimTime/dimLength, scalar(1))
    );

    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Iterate the adjoint to the eikonal equation
    for (label iter = 0; iter < nEikonalIters_; ++iter)
    {
        read();

        Info<< "Adjoint Eikonal Iteration : " << iter << endl;

        fvScalarMatrix daEqn
        (
            2*fvm::div(-yPhi, da_)
          + fvm::SuSp(-epsilon_*fvc::laplacian(d), da_)
          - epsilon_*fvm::laplacian(d, da_)
          + source_
         ==
            fvOptions(scaleDims, da_)
        );

        daEqn.relax();
        fvOptions.constrain(daEqn);
        scalar residual = daEqn.solve().initialResidual();
        fvOptions.correct(da_);

        Info<< "Max da " << gMax(mag(da_)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance_)
        {
            Info<< "\n***Reached adjoint eikonal convergence limit, iteration "
                << iter << "***\n\n";
            break;
        }
    }
    if (debug)
    {
        da_.write();
    }
}


void adjointEikonalSolver::reset()
{
    source_ == dimensionedScalar(source_.dimensions(), Zero);
    distanceSensPtr_() = vector::zero;
}


boundaryVectorField& adjointEikonalSolver::distanceSensitivities()
{
    Info<< "Calculating distance sensitivities " << endl;

    boundaryVectorField& distanceSens = distanceSensPtr_();

    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();

    for (const label patchi : sensitivityPatchIDs_)
    {
        vectorField nf(mesh_.boundary()[patchi].nf());

        // No surface area included. Will be done by the actual sensitivity tool
        distanceSens[patchi] =
           -2.*da_.boundaryField()[patchi]
           *d.boundaryField()[patchi].snGrad()
           *d.boundaryField()[patchi].snGrad()*nf;
    }
    return distanceSens;
}


tmp<volTensorField> adjointEikonalSolver::getFISensitivityTerm() const
{
    Info<< "Calculating distance sensitivities " << endl;

    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();
    const volVectorField gradD(fvc::grad(d));

    auto gradDDa
    (
        tmp<volVectorField>::New
        (
            IOobject
            (
                "gradDDa",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector(d.dimensions()*da_.dimensions()/dimLength, Zero),
            patchDistMethod::patchTypes<vector>(mesh_, wallPatchIDs_)
        )
    );
    gradDDa.ref() = fvc::grad(d*da_);

    auto tdistanceSens
    (
        tmp<volTensorField>::New
        (
            IOobject
            (
                "distanceSensFI",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedTensor(da_.dimensions(), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    volTensorField& distanceSens = tdistanceSens.ref();

    distanceSens =
        - 2.*da_*gradD*gradD
        - epsilon_*gradDDa*gradD
        // grad(gradD) is symmetric theoretically but not numerically when
        // computed with the Gauss divergence theorem. The following maintains
        // exactly the same behaviour as the one before the re-factoring of
        // sensitivities but avoid calling the tranpose operator.
        + epsilon_*da_*d*fvc::div(fvc::interpolate(gradD)*mesh_.Sf());
    distanceSens.correctBoundaryConditions();

    return tdistanceSens;
}


tmp<scalarField> adjointEikonalSolver::topologySensitivities
(
    const word& designVarsName
) const
{
    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();

    auto tres(tmp<scalarField>::New(d.primitiveField().size(), Zero));
    scalarField dSens(d.primitiveField()*da_.primitiveField());

    fv::options& fvOptions(fv::options::New(this->mesh_));
    sensitivityTopO::postProcessSens
    (
        tres.ref(), dSens, fvOptions, d.name(), designVarsName
    );

    return tres;
}


const volScalarField& adjointEikonalSolver::da()
{
    return da_;
}


tmp<volVectorField> adjointEikonalSolver::gradEikonal()
{
    const tmp<volScalarField> td(adjointSolver_.yWall());
    const volScalarField& d = td();

    volVectorField gradD(fvc::grad(d));
    return tmp<volVectorField>::New("gradEikonal", 2*gradD & fvc::grad(gradD));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
