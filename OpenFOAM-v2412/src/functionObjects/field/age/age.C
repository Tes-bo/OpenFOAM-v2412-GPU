/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "age.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvOptions.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "turbulenceModel.H"
#include "inletOutletFvPatchField.H"
#include "wallFvPatch.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(age, 0);
    addToRunTimeSelectionTable(functionObject, age, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::age::patchTypes() const
{
    wordList result
    (
        mesh_.boundary().size(),
        inletOutletFvPatchField<scalar>::typeName
    );

    forAll(mesh_.boundary(), patchi)
    {
        if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
        {
            result[patchi] = fvPatchFieldBase::zeroGradientType();
        }
    }

    return result;
}


bool Foam::functionObjects::age::converged
(
    const int nCorr,
    const scalar initialResidual
) const
{
    if (initialResidual < tolerance_)
    {
        Info<< "Field " << typeName
            << " converged in " << nCorr << " correctors"
            << nl << endl;

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::age::age
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::age::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        phiName_ = dict.getOrDefault<word>("phi", "phi");
        rhoName_ = dict.getOrDefault<word>("rho", "rho");
        schemesField_ = dict.getOrDefault<word>("schemesField", typeName);
        tolerance_ = dict.getOrDefault<scalar>("tolerance", 1e-5);
        nCorr_ = dict.getOrDefault<int>("nCorr", 5);
        diffusion_ = dict.getOrDefault<bool>("diffusion", false);

        return true;
    }

    return false;
}


bool Foam::functionObjects::age::execute()
{
    auto tage = tmp<volScalarField>::New
    (
        IOobject
        (
            typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimTime, Zero),
        patchTypes()
    );
    volScalarField& age = tage.ref();

    const word divScheme("div(phi," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0;
    mesh_.relaxEquation(schemesField_, relaxCoeff);

    Foam::fv::options& fvOptions(Foam::fv::options::New(mesh_));


    // This only works because the null constructed inletValue for an
    // inletOutletFvPatchField is zero. If we needed any other value we would
    // have to loop over the inletOutlet patches and explicitly set the
    // inletValues. We would need to change the interface of inletOutlet in
    // order to do this.

    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimMass/dimTime)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>(rhoName_);

        tmp<volScalarField> tmuEff;
        word laplacianScheme;

        if (diffusion_)
        {
            tmuEff =
                mesh_.lookupObject<compressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                ).muEff();

            laplacianScheme =
                "laplacian(" + tmuEff().name() + ',' + schemesField_ + ")";
        }

        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix ageEqn
            (
                fvm::div(phi, age, divScheme) == rho //+ fvOptions(rho, age)
            );

            if (diffusion_)
            {
                ageEqn -= fvm::laplacian(tmuEff(), age, laplacianScheme);
            }

            ageEqn.relax(relaxCoeff);

            fvOptions.constrain(ageEqn);

            if (converged(i, ageEqn.solve().initialResidual()))
            {
                break;
            };

            fvOptions.correct(age);
        }
    }
    else
    {
        tmp<volScalarField> tnuEff;
        word laplacianScheme;

        if (diffusion_)
        {
            tnuEff =
                mesh_.lookupObject<incompressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                ).nuEff();

            laplacianScheme =
                "laplacian(" + tnuEff().name() + ',' + schemesField_ + ")";
        }

        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix ageEqn
            (
                fvm::div(phi, age, divScheme)
             == dimensionedScalar(1) + fvOptions(age)
            );

            if (diffusion_)
            {
                ageEqn -= fvm::laplacian(tnuEff(), age, laplacianScheme);
            }

            ageEqn.relax(relaxCoeff);

            fvOptions.constrain(ageEqn);

            if (converged(i, ageEqn.solve().initialResidual()))
            {
                break;
            }

            fvOptions.correct(age);
        }
    }

    Info<< "Min/max age:"
        << min(age).value() << ' '
        << max(age).value()
        << endl;

    // Workaround
    word fieldName = typeName;

    return store(fieldName, tage);
}


bool Foam::functionObjects::age::write()
{
    return writeObject(typeName);
}


// ************************************************************************* //
