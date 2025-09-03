/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "advectionDiffusionPatchDistMethod.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(advectionDiffusion, 0);
    addToRunTimeSelectionTable(patchDistMethod, advectionDiffusion, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::advectionDiffusion::advectionDiffusion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    //- We do not want to recurse into 'advectionDiffusion' again so
    //  make sure we pick up 'method' always from the subdictionary.
    //coeffs_(dict.optionalSubDict(type() + "Coeffs")),
    coeffs_(dict.subDict(type() + "Coeffs")),
    pdmPredictor_
    (
        patchDistMethod::New
        (
            coeffs_,
            mesh,
            patchIDs
        )
    ),
    epsilon_(coeffs_.getOrDefault<scalar>("epsilon", 0.1)),
    tolerance_(coeffs_.getOrDefault<scalar>("tolerance", 1e-3)),
    maxIter_(coeffs_.getOrDefault<int>("maxIter", 10)),
    predicted_(false),
    checkAndWriteMesh_(coeffs_.getOrDefault("checkAndWriteMesh", false))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::advectionDiffusion::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::advectionDiffusion::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    if (!predicted_)
    {
        pdmPredictor_->correct(y);
        predicted_ = true;
    }

    // If the mesh has become invalid, trying to solve the eikonal equation
    // might lead to crashing and we wont get a chance to see the failed mesh.
    // Write the mesh points here
    if (checkAndWriteMesh_)
    {
        DebugInfo
            << "Checking mesh in advectionDiffusion " << endl;
        if (mesh_.checkMesh(true))
        {
            pointIOField points
            (
                IOobject
                (
                   "points",
                    mesh_.pointsInstance(),
                    polyMesh::meshSubDir,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_.points()
            );
            points.write();
        }
    }

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
        patchTypes<vector>(mesh_, patchIDs_)
    );

    const fvPatchList& patches = mesh_.boundary();
    volVectorField::Boundary& nybf = ny.boundaryFieldRef();

    for (const label patchi : patchIDs_)
    {
        nybf[patchi] == -patches[patchi].nf();
    }

    // Scaling dimensions, to be used with fvOptions
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
    int iter = 0;
    scalar initialResidual = 0;

    do
    {
        ny = fvc::grad(y);
        ny /= (mag(ny) + SMALL);

        surfaceVectorField nf(fvc::interpolate(ny));
        nf /= (mag(nf) + SMALL);

        surfaceScalarField yPhi("yPhi", mesh_.Sf() & nf);

        fvScalarMatrix yEqn
        (
            fvm::div(yPhi, y)
          + fvm::SuSp(-fvc::div(yPhi), y)
          - epsilon_*y*fvm::laplacian(y)
         ==
            dimensionedScalar("1", dimless, 1.0)
          + fvOptions(scaleDims, y)
        );

        yEqn.relax();
        fvOptions.constrain(yEqn);
        initialResidual = yEqn.solve().initialResidual();
        fvOptions.correct(y);

    } while (initialResidual > tolerance_ && ++iter < maxIter_);

    // Need to stabilise the y for overset meshes since the holed cells
    // keep the initial value (0.0) so the gradient of that will be
    // zero as well. Turbulence models do not like zero wall distance.
    y.clamp_min(SMALL);

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n = -ny;
    }

    return true;
}


// ************************************************************************* //
