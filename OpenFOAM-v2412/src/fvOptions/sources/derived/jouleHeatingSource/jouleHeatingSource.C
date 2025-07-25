/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "jouleHeatingSource.H"
#include "fvMatrices.H"
#include "fvmLaplacian.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(jouleHeatingSource, 0);
    addToRunTimeSelectionTable(option, jouleHeatingSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::fv::jouleHeatingSource::transformSigma
(
    const volVectorField& sigmaLocal
) const
{
    auto tsigma = volSymmTensorField::New
    (
        IOobject::scopedName(typeName, "sigma"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedSymmTensor(sigmaLocal.dimensions(), Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    auto& sigma = tsigma.ref();

    // This check should be unnecessary
    if (!csysPtr_)
    {
        FatalErrorInFunction
            << "Coordinate system undefined"
            << abort(FatalError);
    }

    const auto& csys = *csysPtr_;

    if (csys.uniform())
    {
        sigma.primitiveFieldRef() =
            csys.transformPrincipal(sigmaLocal);
    }
    else
    {
        sigma.primitiveFieldRef() =
            csys.transformPrincipal(mesh_.cellCentres(), sigmaLocal);
    }

    sigma.correctBoundaryConditions();

    return tsigma;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::jouleHeatingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    TName_("T"),
    V_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "V"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh
    ),
    anisotropicElectricalConductivity_(false),
    scalarSigmaVsTPtr_(nullptr),
    vectorSigmaVsTPtr_(nullptr),
    csysPtr_(nullptr),
    curTimeIndex_(-1)
{
    // Set the field name to that of the energy
    // field from which the temperature is obtained

    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::fv::jouleHeatingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo<< name() << ": applying source to " << eqn.psi().name() << endl;

    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        if (anisotropicElectricalConductivity_)
        {
            // Update sigma as a function of T if required
            const volVectorField& sigmaLocal = updateSigma(vectorSigmaVsTPtr_);

            tmp<volSymmTensorField> sigma = transformSigma(sigmaLocal);

            // Solve the electrical potential equation
            fvScalarMatrix VEqn(fvm::laplacian(sigma, V_));
            VEqn.relax();
            VEqn.solve();
        }
        else
        {
            // Update sigma as a function of T if required
            const volScalarField& sigma = updateSigma(scalarSigmaVsTPtr_);

            // Solve the electrical potential equation
            fvScalarMatrix VEqn(fvm::laplacian(sigma, V_));
            VEqn.relax();
            VEqn.solve();
        }

        curTimeIndex_ = mesh_.time().timeIndex();
    }

    // Add the Joule heating contribution
    const word sigmaName
    (
        IOobject::scopedName(typeName, "sigma")
    );

    const volVectorField gradV(fvc::grad(V_));
    if (anisotropicElectricalConductivity_)
    {
        const auto& sigmaLocal = mesh_.lookupObject<volVectorField>(sigmaName);

        tmp<volSymmTensorField> sigma = transformSigma(sigmaLocal);

        eqn += (sigma & gradV) & gradV;
    }
    else
    {
        const auto& sigma = mesh_.lookupObject<volScalarField>(sigmaName);

        eqn += (sigma*gradV) & gradV;
    }
}


bool Foam::fv::jouleHeatingSource::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        coeffs_.readIfPresent("T", TName_);

        anisotropicElectricalConductivity_ =
            coeffs_.get<bool>("anisotropicElectricalConductivity");

        if (anisotropicElectricalConductivity_)
        {
            Info<< "    Using vector electrical conductivity" << endl;

            initialiseSigma(coeffs_, vectorSigmaVsTPtr_);

            csysPtr_ =
                coordinateSystem::New
                (
                    mesh_,
                    coeffs_,
                    coordinateSystem::typeName
                );
        }
        else
        {
            Info<< "    Using scalar electrical conductivity" << endl;

            initialiseSigma(coeffs_, scalarSigmaVsTPtr_);

            csysPtr_.clear();  // Do not need coordinate system
        }

        return true;
    }

    return false;
}

// ************************************************************************* //
