/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "scalarTransport.H"
#include "CMULES.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(scalarTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        scalarTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::scalarTransport::transportedField()
{
    if (!foundObject<volScalarField>(fieldName_))
    {
        auto tfldPtr = tmp<volScalarField>::New
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_
        );
        store(fieldName_, tfldPtr);

        if (phaseName_ != "none")
        {
            mesh_.setFluxRequired(fieldName_);
        }
    }

    return lookupObjectRef<volScalarField>(fieldName_);
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::scalarTransport::D
(
    const volScalarField& s,
    const surfaceScalarField& phi
) const
{
    const word Dname("D" + s.name());

    if (constantD_)
    {
        return volScalarField::New
        (
            Dname,
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(Dname, phi.dimensions()/dimLength, D_)
        );
    }

    if (nutName_ != "none")
    {
        const volScalarField& nutMean =
            mesh_.lookupObject<volScalarField>(nutName_);

        return tmp<volScalarField>::New(Dname, nutMean);
    }

    // Incompressible
    {
        const auto* turb =
            findObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turb)
        {
            return volScalarField::New
            (
                Dname,
                IOobject::NO_REGISTER,
                alphaD_ * turb->nu() + alphaDt_ * turb->nut()
            );
        }
    }

    // Compressible
    {
        const auto* turb =
            findObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turb)
        {
            return volScalarField::New
            (
                Dname,
                IOobject::NO_REGISTER,
                alphaD_ * turb->mu() + alphaDt_ * turb->mut()
            );
        }
    }


    return volScalarField::New
    (
        Dname,
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(phi.dimensions()/dimLength, Zero)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransport::scalarTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fvOptions_(mesh_),
    fieldName_(dict.getOrDefault<word>("field", "s")),
    schemesField_("unknown-schemesField"),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    nutName_(dict.getOrDefault<word>("nut", "none")),
    phaseName_(dict.getOrDefault<word>("phase", "none")),
    phasePhiCompressedName_
    (
        dict.getOrDefault<word>("phasePhiCompressed", "alphaPhiUn")
    ),
    D_(0),
    alphaD_(1),
    alphaDt_(1),
    tol_(1),
    nCorr_(0),
    resetOnStartUp_(false),
    constantD_(false),
    bounded01_(dict.getOrDefault("bounded01", true))
{
    read(dict);

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& s = transportedField();

    if (resetOnStartUp_)
    {
        s == dimensionedScalar(dimless, Zero);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::scalarTransport::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    dict.readIfPresent("phi", phiName_);
    dict.readIfPresent("rho", rhoName_);
    dict.readIfPresent("nut", nutName_);
    dict.readIfPresent("phase", phaseName_);
    dict.readIfPresent("phasePhiCompressed", phasePhiCompressedName_);

    schemesField_ = dict.getOrDefault("schemesField", fieldName_);

    dict.readIfPresent("alphaD", alphaD_);
    dict.readIfPresent("alphaDt", alphaDt_);
    dict.readIfPresent("tolerance", tol_);
    dict.readIfPresent("nCorr", nCorr_);
    dict.readIfPresent("resetOnStartUp", resetOnStartUp_);
    constantD_ = dict.readIfPresent("D", D_);
    dict.readIfPresent("bounded01", bounded01_);

    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    return true;
}


bool Foam::functionObjects::scalarTransport::execute()
{
    volScalarField& s = transportedField();

    Log << type() << " execute: " << s.name() << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Calculate the diffusivity
    volScalarField D(this->D(s, phi));

    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0;
    mesh_.relaxEquation(schemesField_, relaxCoeff);

    // Convergence monitor parameters
    bool converged = false;
    int iter = 0;

    // Two phase scalar transport
    if (phaseName_ != "none")
    {
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(phaseName_);

        const surfaceScalarField& limitedPhiAlpha =
            mesh_.lookupObject<surfaceScalarField>(phasePhiCompressedName_);

        D *= pos(alpha - 0.99);

        // Reset D dimensions consistent with limitedPhiAlpha
        D.dimensions().reset(limitedPhiAlpha.dimensions()/dimLength);

        // Solve
        tmp<surfaceScalarField> tTPhiUD;
        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s)
              + fvm::div(limitedPhiAlpha, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
              ==
                alpha*fvOptions_(s)
            );

            sEqn.relax(relaxCoeff);
            fvOptions_.constrain(sEqn);

            ++iter;
            converged = (sEqn.solve(schemesField_).initialResidual() < tol_);

            tTPhiUD = sEqn.flux();

            if (converged) break;
        }

        if (bounded01_)
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                s,
                phi,
                tTPhiUD.ref(),
                oneField(),
                zeroField()
            );
        }
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho = lookupObject<volScalarField>(rhoName_);

        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s)
              + fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions_(rho, s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            ++iter;
            converged = (sEqn.solve(schemesField_).initialResidual() < tol_);
            if (converged) break;
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s)
              + fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions_(s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            ++iter;
            converged = (sEqn.solve(schemesField_).initialResidual() < tol_);
            if (converged) break;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    if (converged)
    {
        Log << type() << ": " << name() << ": "
            << s.name() << " is converged." << nl
            << tab << "initial-residual tolerance: " << tol_ << nl
            << tab << "outer iteration: " << iter << nl;
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::scalarTransport::write()
{
    Log << type() << " write: " << name() << nl
        << tab << fieldName_ << nl
        << endl;

    volScalarField& s = transportedField();

    s.write();

    return true;
}


// ************************************************************************* //
