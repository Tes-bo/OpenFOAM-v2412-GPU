/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "phaseModel.H"
#include "diameterModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const word& phaseName,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    name_(phaseName),
    phaseDict_(phaseDict),
    nu_
    (
        "nu",
        dimViscosity,
        phaseDict_
    ),
    kappa_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        phaseDict_
    ),
    Cp_
    (
        "Cp",
        dimSpecificHeatCapacity,
        phaseDict_
    ),
    rho_
    (
        "rho",
        dimDensity,
        phaseDict_
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    DDtU_
    (
        IOobject
        (
            IOobject::groupName("DDtU", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    )
{
    alphaPhi_.setOriented();

    IOobject io
    (
        IOobject::groupName("phi", name_),
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );

    if (io.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << io.name() << endl;

        phiPtr_.reset(new surfaceScalarField(io, mesh));
    }
    else
    {
        Info<< "Calculating face flux field " << io.name() << endl;

        io.readOpt(IOobject::NO_READ);

        wordList phiTypes
        (
            U_.boundaryField().size(),
            fvsPatchFieldBase::calculatedType()
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                io,
                fvc::flux(U_),
                phiTypes
            )
        );
    }

    dPtr_ = multiphaseEuler::diameterModel::New
    (
        phaseDict_,
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return nullptr;
}



void Foam::phaseModel::correct()
{
    //nuModel_->correct();
}


bool Foam::phaseModel::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    //if (nuModel_->read(phaseDict_))
    {
        phaseDict_.readEntry("nu", nu_.value());
        phaseDict_.readEntry("kappa", kappa_.value());
        phaseDict_.readEntry("Cp", Cp_.value());
        phaseDict_.readEntry("rho", rho_.value());

        return true;
    }

    // return false;
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const tmp<surfaceScalarField> tphi(phi());
    const auto& phiBf = tphi().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}


// ************************************************************************* //
