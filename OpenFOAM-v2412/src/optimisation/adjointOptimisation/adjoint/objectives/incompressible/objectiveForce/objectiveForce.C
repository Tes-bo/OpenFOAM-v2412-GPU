/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023, 2022 PCOpt/NTUA
    Copyright (C) 2013-2023, 2022 FOSS GP
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "objectiveForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveForce, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveForce::objectiveForce
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    forcePatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    ),
    forceDirection_(dict.get<vector>("direction")),
    Aref_(dict.get<scalar>("Aref")),
    rhoInf_(dict.get<scalar>("rhoInf")),
    UInf_(dict.get<scalar>("UInf"))
{
    // Sanity check and print info
    if (forcePatches_.empty())
    {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }
    if (debug)
    {
        Info<< "Minimizing " << type() << " in patches:" << endl;
        for (const label patchI : forcePatches_)
        {
            Info<< "\t " << mesh_.boundary()[patchI].name() << endl;
        }
    }

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdnutPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    bdJdGradUPtr_.reset(createZeroBoundaryPtr<tensor>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveForce::J()
{
    vector pressureForce(Zero);
    vector viscousForce(Zero);
    vector cumulativeForce(Zero);


    const volScalarField& p = vars_.pInst();
    const autoPtr<incompressible::turbulenceModel>&
       turbulence = vars_.turbulence();

    volSymmTensorField devReff(turbulence->devReff());

    for (const label patchI : forcePatches_)
    {
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        pressureForce += gSum(Sf*p.boundaryField()[patchI]);
        viscousForce += gSum(devReff.boundaryField()[patchI] & Sf);
    }

    cumulativeForce = pressureForce + viscousForce;

    scalar force = cumulativeForce & forceDirection_;

    // Intentionally not using denom - derived might implement virtual denom()
    // function differently
    scalar Cforce = force/(0.5*UInf_*UInf_*Aref_);

    DebugInfo
        << "Force|Coeff " << force << "|" << Cforce << endl;

    J_ = Cforce;

    return Cforce;
}


void objectiveForce::update_boundarydJdp()
{
    for (const label patchI : forcePatches_)
    {
        bdJdpPtr_()[patchI] = forceDirection_/denom();
    }
}


void objectiveForce::update_dSdbMultiplier()
{
    // Compute contributions with mean fields, if present
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();
    const autoPtr<incompressible::RASModelVariables>& turbVars =
        vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

    tmp<volSymmTensorField> tdevReff = turbVars->devReff(lamTransp, U);
    const volSymmTensorField& devReff = tdevReff();

    for (const label patchI : forcePatches_)
    {
        bdSdbMultPtr_()[patchI] =
        (
            (
                forceDirection_& devReff.boundaryField()[patchI]
            )
          + (forceDirection_)*p.boundaryField()[patchI]
        )
       /denom();
    }
}


void objectiveForce::update_dxdbMultiplier()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    const autoPtr<incompressible::RASModelVariables>&
        turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();

    // We only need to modify the boundaryField of gradU locally.
    // If grad(U) is cached then
    // a. The .ref() call fails since the tmp is initialised from a
    //    const ref
    // b. we would be changing grad(U) for all other places in the code
    //    that need it
    // So, always allocate new memory and avoid registering the new field
    tmp<volTensorField> tgradU =
        volTensorField::New("gradULocal", fvc::grad(U));
    volTensorField& gradU = tgradU.ref();
    volTensorField::Boundary& gradUbf = gradU.boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = patch.nf();
            gradUbf[patchI] = tnf*U.boundaryField()[patchI].snGrad();
        }
    }

    // Term coming from gradp
    tmp<volVectorField> tgradp(fvc::grad(p));
    const volVectorField& gradp = tgradp.cref();
    for (const label patchI : forcePatches_)
    {
        bdxdbMultPtr_()[patchI] =
            (forceDirection_ & mesh_.boundary()[patchI].nf())
           *gradp.boundaryField()[patchI]/denom();
    }
    tgradp.clear();

    // Term coming from stresses
    tmp<volScalarField> tnuEff = lamTransp.nu() + turbVars->nut();
    tmp<volSymmTensorField> tstress = tnuEff*twoSymm(tgradU);
    const volSymmTensorField& stress = tstress.cref();
    autoPtr<volVectorField> ptemp
        (Foam::createZeroFieldPtr<vector>( mesh_, "temp", sqr(dimVelocity)));
    volVectorField& temp = ptemp.ref();

    for (label idir = 0; idir < pTraits<vector>::nComponents; ++idir)
    {
        unzipRow(stress, idir, temp);
        volTensorField gradStressDir(fvc::grad(temp));
        for (const label patchI : forcePatches_)
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            tmp<vectorField> tnf = patch.nf();
            bdxdbMultPtr_()[patchI] -=
                forceDirection_.component(idir)
               *(gradStressDir.boundaryField()[patchI] & tnf)/denom();
        }
    }
}


void objectiveForce::update_boundarydJdnut()
{
    const volVectorField& U = vars_.U();
    volSymmTensorField devGradU(devTwoSymm(fvc::grad(U)));

    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        bdJdnutPtr_()[patchI] =
          - ((devGradU.boundaryField()[patchI] & forceDirection_) & tnf)
           /denom();
    }
}


void objectiveForce::update_boundarydJdGradU()
{
    const autoPtr<incompressible::RASModelVariables>& turbVars =
        vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nut());
    for (const label patchI : forcePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const vectorField& Sf = patch.Sf();
        bdJdGradUPtr_()[patchI] =
          - nuEff.boundaryField()[patchI]
           *dev(forceDirection_*Sf + Sf*forceDirection_)/denom();
    }
}


scalar objectiveForce::denom() const
{
    return 0.5*UInf_*UInf_*Aref_;
}


const vector& objectiveForce::forceDirection() const
{
    return forceDirection_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
