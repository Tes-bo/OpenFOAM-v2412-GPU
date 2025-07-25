/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "electricPotential.H"
#include "fvc.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(electricPotential, 0);
    addToRunTimeSelectionTable(functionObject, electricPotential, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::electricPotential::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = mesh_.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::electricPotential::sigma() const
{
    const IOobject sigmaIO
    (
        mesh_.thisDb().newIOobject(IOobject::scopedName(typeName, "sigma"))
    );

    if (phases_.size())
    {
        tmp<volScalarField> tsigma = phases_[0]*sigmas_[0];

        for (label i = 1; i < phases_.size(); ++i)
        {
            tsigma.ref() += phases_[i]*sigmas_[i];
        }

        return tmp<volScalarField>::New
        (
            sigmaIO,
            tsigma,
            fvPatchFieldBase::calculatedType()
        );
    }

    return tmp<volScalarField>::New
    (
        sigmaIO,
        mesh_,
        sigma_,
        fvPatchFieldBase::calculatedType()
    );
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::electricPotential::epsilonm() const
{
    // Vacuum permittivity (aka the electric constant) [A^2 s^4/(kg m^3)]
    const dimensionedScalar epsilon0
    (
        sqr(dimCurrent)*pow4(dimTime)/(dimMass*pow3(dimLength)),
        8.8541878128e-12    // CODATA value
    );

    const IOobject epsilonrIO
    (
        mesh_.thisDb().newIOobject(IOobject::scopedName(typeName, "epsilonr"))
    );

    if (phases_.size())
    {
        tmp<volScalarField> tepsilonr = phases_[0]*epsilonrs_[0];

        for (label i = 1; i < phases_.size(); ++i)
        {
            tepsilonr.ref() += phases_[i]*epsilonrs_[i];
        }

        return tmp<volScalarField>::New
        (
            epsilonrIO,
            epsilon0*tepsilonr,
            fvPatchFieldBase::calculatedType()
        );
    }

    return tmp<volScalarField>::New
    (
        epsilonrIO,
        mesh_,
        epsilon0*epsilonr_,
        fvPatchFieldBase::calculatedType()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::electricPotential::electricPotential
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phasesDict_(dict.subOrEmptyDict("phases")),
    phaseNames_(),
    phases_(),
    sigmas_(),
    sigma_
    (
        dimensionedScalar
        (
            sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
            dict.getCheckOrDefault<scalar>
            (
                "sigma",
                scalar(1),
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    epsilonrs_(),
    epsilonr_
    (
        dimensionedScalar
        (
            dimless,
            dict.getCheckOrDefault<scalar>
            (
                "epsilonr",
                scalar(1),
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    Vname_
    (
        dict.getOrDefault<word>
        (
            "V",
            IOobject::scopedName(typeName, "V")
        )
    ),
    Ename_
    (
        dict.getOrDefault<word>
        (
            "E",
            IOobject::scopedName(typeName, "E")
        )
    ),
    fvOptions_(mesh_),
    tol_(1),
    nCorr_(1),
    writeDerivedFields_(false),
    electricField_(false)
{
    read(dict);

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& eV = getOrReadField(Vname_);
    eV.correctBoundaryConditions();

    if (electricField_)
    {
        auto* ptr = getObjectPtr<volVectorField>(Ename_);

        if (!ptr)
        {
            ptr = new volVectorField
            (
                IOobject
                (
                    Ename_,
                    mesh_.time().timeName(),
                    mesh_.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                -fvc::grad(eV)
            );
            regIOobject::store(ptr);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::electricPotential::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    Log << type() << " read: " << name() << endl;

    dict.readIfPresent("sigma", sigma_);
    dict.readIfPresent("epsilonr", epsilonr_);
    dict.readIfPresent("nCorr", nCorr_);
    dict.readIfPresent("tolerance", tol_);
    dict.readIfPresent("writeDerivedFields", writeDerivedFields_);
    dict.readIfPresent("electricField", electricField_);

    // If flow is multiphase
    if (!phasesDict_.empty())
    {
        phaseNames_.setSize(phasesDict_.size());
        phases_.setSize(phasesDict_.size());
        sigmas_.setSize(phasesDict_.size());

        if (writeDerivedFields_)
        {
            epsilonrs_.setSize(phasesDict_.size());
        }

        label phasei = 0;
        for (const entry& dEntry : phasesDict_)
        {
            const word& key = dEntry.keyword();

            if (!dEntry.isDict())
            {
                FatalIOErrorInFunction(phasesDict_)
                    << "Entry " << key << " is not a dictionary" << nl
                    << exit(FatalIOError);
            }

            const dictionary& subDict = dEntry.dict();

            phaseNames_[phasei] = key;

            sigmas_.set
            (
                phasei,
                new dimensionedScalar
                (
                    sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
                    subDict.getCheck<scalar>
                    (
                        "sigma",
                        scalarMinMax::ge(SMALL)
                    )
                )
            );

            if (writeDerivedFields_)
            {
                epsilonrs_.set
                (
                    phasei,
                    new dimensionedScalar
                    (
                        dimless,
                        subDict.getCheck<scalar>
                        (
                            "epsilonr",
                            scalarMinMax::ge(SMALL)
                        )
                    )
                );
            }

            ++phasei;
        }

        forAll(phaseNames_, i)
        {
            phases_.set
            (
                i,
                mesh_.getObjectPtr<volScalarField>(phaseNames_[i])
            );
        }
    }

    if (const dictionary* dictptr = dict.findDict("fvOptions"))
    {
        fvOptions_.reset(*dictptr);
    }

    return true;
}


bool Foam::functionObjects::electricPotential::execute()
{
    Log << type() << " execute: " << name() << endl;

    // Convergence monitor parameters
    bool converged = false;
    label iter = 0;

    tmp<volScalarField> tsigma = this->sigma();
    const auto& sigma = tsigma();

    volScalarField& eV = getOrReadField(Vname_);

    for (int i = 1; i <= nCorr_; ++i)
    {
        fvScalarMatrix eVEqn
        (
          - fvm::laplacian(sigma, eV)
        );

        eVEqn.relax();

        fvOptions_.constrain(eVEqn);

        ++iter;
        converged = (eVEqn.solve().initialResidual() < tol_);
        if (converged) break;
    }

    if (electricField_)
    {
        auto& E = lookupObjectRef<volVectorField>(Ename_);
        E == -fvc::grad(eV);
    }

    if (converged)
    {
        Log << type() << ": " << name() << ": "
            << eV.name() << " is converged." << nl
            << tab << "initial-residual tolerance: " << tol_ << nl
            << tab << "outer iteration: " << iter << nl;
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::electricPotential::write()
{
    Log << type() << " write: " << name() << nl
        << tab << Vname_
        << endl;

    volScalarField& eV = getOrReadField(Vname_);

    if (electricField_)
    {
        const auto& E = lookupObject<volVectorField>(Ename_);

        Log << tab << E.name() << endl;

        E.write();
    }

    if (writeDerivedFields_)
    {
        // Write the current density field
        tmp<volScalarField> tsigma = this->sigma();

        auto eJ = volVectorField::New
        (
            IOobject::scopedName(typeName, "J"),
            IOobject::NO_REGISTER,
            -tsigma*fvc::grad(eV),
            fvPatchFieldBase::calculatedType()
        );

        Log << tab << eJ().name() << endl;

        eJ->write();


        // Write the free-charge density field
        tmp<volScalarField> tepsilonm = this->epsilonm();

        auto erho = volScalarField::New
        (
            IOobject::scopedName(typeName, "rho"),
            IOobject::NO_REGISTER,
            fvc::div(tepsilonm*(-fvc::grad(eV))),
            fvPatchFieldBase::calculatedType()
        );

        Log << tab << erho().name() << endl;

        erho->write();
    }

    eV.write();

    return true;
}


// ************************************************************************* //
