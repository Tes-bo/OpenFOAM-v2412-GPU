/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

Application
    setTurbulenceFields

Group
    grpPreProcessingUtilities

Description
    Initialises turbulence fields according to
    various empirical governing equations.

    References:
    \verbatim
        Initialisation method (tag:M):
            Manceau, R. (n.d.).
            A two-step automatic initialization
            procedure for RANS computations.
            (Unpublished).

        Previous-version of the initialisation model (tag:LM):
            Lardeau, S., & Manceau, R. (2014).
            Computations of complex flow configurations using
            a modified elliptic-blending Reynolds-stress model.
            10th International ERCOFTAC Symposium on Engineering
            Turbulence Modelling and Measurements. Marbella, Spain.
            https://hal.archives-ouvertes.fr/hal-01051799
    \endverbatim

Usage
    Minimal example by using \c system/setTurbulenceFieldsDict:
    \verbatim
        // Mandatory entries
        uRef            <scalar>;

        // Optional entries
        initialiseU     <bool>;
        initialiseEpsilon <bool>;
        initialiseK     <bool>;
        initialiseOmega <bool>;
        initialiseR     <bool>;
        writeF          <bool>;

        kappa           <scalar>;
        Cmu             <scalar>;
        dPlusRef        <scalar>;

        f               <word>;
        U               <word>;
        epsilon         <word>;
        k               <word>;
        omega           <word>;
        R               <word>;
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type | Reqd | Deflt
      uRef      | Reference speed              | scalar | yes  | -
      initialiseU | Flag to initialise U         | bool | no   | false
      initialiseEpsilon | Flag to initialise epsilon    | bool | no   | false
      initialiseK | Flag to initialise k         | bool | no   | false
      initialiseOmega | Flag to initialise omega | bool | no   | false
      initialiseR | Flag to initialise R         | bool | no   | false
      writeF | Flag to write elliptic-blending field, f | bool | no   | false
      kappa  | von Karman constant             | scalar | no   | 0.41
      Cmu    | Empirical constant              | scalar | no   | 0.09
      dPlusRef | Reference dPlus               | scalar | no   | 15
      f      | Name of operand f field         | word   | no   | f
      U      | Name of operand U field         | word   | no   | U
      epsilon  | Name of operand epsilon field | word   | no   | epsilon
      k      | Name of operand k field         | word   | no   | k
      omega      | Name of operand omega field | word   | no   | omega
      R      | Name of operand R field         | word   | no   | R
    \endtable

Note
  - Time that the utility applies to is determined by the
    \c controlDict.startFrom and \c controlDict.startTime entries.
  - The utility modifies near-wall fields, hence
    can be more effective for low-Re mesh cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallFvPatch.H"
#include "processorFvPatch.H"
#include "fixedValueFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void InfoField(const word& fldName)
{
    Info<< "Writing field: " << fldName << nl << endl;
}


template<class Type>
void correctProcessorPatches(GeometricField<Type, fvPatchField, volMesh>& fld)
{
    if (UPstream::parRun())
    {
        fld.boundaryFieldRef().template evaluateCoupled<processorFvPatch>();
    }
}


IOobject createIOobject
(
    const fvMesh& mesh,
    const word& name,
    IOobject::readOption rOpt = IOobject::READ_IF_PRESENT
)
{
    return
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            rOpt,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );
}


tmp<volScalarField> nu
(
    const fvMesh& mesh,
    const volVectorField& U
)
{
    const Time& runTime = mesh.time();

    if
    (
        IOobject
        (
            basicThermo::dictName,
            runTime.constant(),
            mesh
        ).typeHeaderOk<IOdictionary>(true)
    )
    {
        // Compressible
        autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
        fluidThermo& thermo = pThermo();
        volScalarField rho(thermo.rho());

        // Create phi
        #include "compressibleCreatePhi.H"

        autoPtr<compressible::turbulenceModel> turbulence
        (
            compressible::turbulenceModel::New
            (
                rho,
                U,
                phi,
                thermo
            )
        );

        turbulence->validate();

        return tmp<volScalarField>::New(turbulence->nu());
    }
    else
    {
        // Incompressible

        // Create phi
        #include "createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
            incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );

        turbulence->validate();

        return tmp<volScalarField>::New(turbulence->nu());
    }
}


// Calculate elliptic blending function
// between near-wall and weakly-inhomogeneous regions
void calcF
(
    const volScalarField& L,
    volScalarField& f
)
{
    tmp<volScalarField> tinvLsqr = scalar(1)/sqr(L);
    const volScalarField& invLsqr = tinvLsqr.cref();

    // (M:Eq. 6)
    tmp<fvScalarMatrix> fEqn
    (
        fvm::Sp(invLsqr, f)
      - fvm::laplacian(f)
      ==
        invLsqr
    );

    tinvLsqr.clear();

    fEqn.ref().relax();
    solve(fEqn);

    // (M:p. 2)
    f.clamp_range(0, scalar(1) - Foam::exp(-scalar(400)/scalar(50)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Sets initial turbulence fields based on"
        " various empirical equations"
    );
    argList::noFunctionObjects();
    argList::addOption("dict", "file", "Alternative setTurbulenceFieldsDict");

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("setTurbulenceFieldsDict");
    #include "setSystemRunTimeDictionaryIO.H"
    Info<< "Reading " << dictIO.name() << nl << endl;
    IOdictionary dict(dictIO);

    #include "createNamedMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOstream::defaultPrecision(15);


    // Dictionary input (M:p. 2)

    const scalar uRef = dict.getCheck<scalar>("uRef", scalarMinMax::ge(SMALL));
    const dimensionedScalar uTau(dimVelocity, 0.05*uRef);
    const scalar kappa =
        dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(SMALL));
    const scalar Cmu =
        dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL));
    const scalar dPlusRef =
        dict.getCheckOrDefault<scalar>("dPlusRef", 15, scalarMinMax::ge(SMALL));
    const word fName = dict.getOrDefault<word>("f", "f");
    const word UName = dict.getOrDefault<word>("U", "U");
    const word epsilonName = dict.getOrDefault<word>("epsilon", "epsilon");
    const word kName = dict.getOrDefault<word>("k", "k");
    const word omegaName = dict.getOrDefault<word>("omega", "omega");
    const word RName = dict.getOrDefault<word>("R", "R");
    const bool initU = dict.getOrDefault<bool>("initialiseU", false);
    const bool initEpsilon =
        dict.getOrDefault<bool>("initialiseEpsilon", false);
    const bool initK = dict.getOrDefault<bool>("initialiseK", false);
    const bool initOmega = dict.getOrDefault<bool>("initialiseOmega", false);
    const bool initR = dict.getOrDefault<bool>("initialiseR", false);
    const bool writeF = dict.getOrDefault<bool>("writeF", false);


    // Start initialising the operand fields

    // Read operand fields

    auto tU = tmp<volVectorField>::New
    (
        createIOobject(mesh, UName, IOobject::MUST_READ),
        mesh
    );

    // Infer the initial BCs from the velocity
    const wordList bcTypes
    (
        tU.cref().boundaryField().size(),
        fixedValueFvPatchScalarField::typeName
    );

    tmp<volScalarField> tepsilon;
    tmp<volScalarField> tk;
    tmp<volScalarField> tomega;
    tmp<volSymmTensorField> tR;

    if (initEpsilon)
    {
        tepsilon = tmp<volScalarField>::New
        (
            createIOobject(mesh, epsilonName),
            mesh,
            dimensionedScalar(sqr(dimLength)/pow3(dimTime), SMALL),
            bcTypes
        );
    }

    if (initK)
    {
        tk = tmp<volScalarField>::New
        (
            createIOobject(mesh, kName),
            mesh,
            dimensionedScalar(sqr(dimLength/dimTime), SMALL),
            bcTypes
        );
    }

    if (initOmega)
    {
        tomega = tmp<volScalarField>::New
        (
            createIOobject(mesh, omegaName),
            mesh,
            dimensionedScalar(dimless/dimTime, SMALL),
            bcTypes
        );
    }

    if (initR)
    {
        tR = tmp<volSymmTensorField>::New
        (
            createIOobject(mesh, RName),
            mesh,
            dimensionedSymmTensor(sqr(dimLength/dimTime), Zero),
            bcTypes
        );
    }


    // Create elliptic blending factor field

    volScalarField f
    (
        IOobject
        (
            fName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, scalar(1)),
        fixedValueFvPatchScalarField::typeName
    );

    for (fvPatchScalarField& pfld : f.boundaryFieldRef())
    {
        if (isA<wallFvPatch>(pfld.patch()))
        {
            pfld == Zero;
        }
    }


    // Create auxillary fields for the initialisation

    tmp<volScalarField> tnu = nu(mesh, tU.cref());
    const volScalarField& nu = tnu.cref();

    // (M:p. 2)
    tmp<volScalarField> tL = scalar(50)*nu/uTau;
    const volScalarField& L = tL.cref();

    calcF(L, f);

    // (M:Eq. 8)
    tmp<volScalarField> td = -tL*Foam::log(scalar(1) - f);
    const volScalarField& d = td.cref();

    // (M:p. 2)
    const volScalarField dPlus(d*uTau/nu);

    // (M:Eq. 11)
    const volScalarField epsilon
    (
        pow4(uTau)/(kappa*nu*max(dPlus, dPlusRef))
    );

    // (M:Eq. 13)
    const volScalarField fk(Foam::exp(-dPlus/scalar(25)));

    // (M:Eq. 12)
    const volScalarField k
    (
        (epsilon*sqr(td)*sqr(fk))/(2*nu)
        + sqr(uTau)*sqr(scalar(1) - fk)/Foam::sqrt(Cmu)
    );


    // Initialise operand fields

    if (initU)
    {
        volVectorField& U = tU.ref();

        // Reichardt’s law (M:Eq. 10)
        const scalar C = 7.8;
        const scalar B1 = 11;
        const scalar B2 = 3;
        const volScalarField fRei
        (
            Foam::log(scalar(1) + kappa*dPlus)/kappa
          + C*
            (
                scalar(1)
              - Foam::exp(-dPlus/B1)
              - dPlus/B1*Foam::exp(-dPlus/B2)
            )
        );

        // (M:Eq. 9)
        const dimensionedScalar maxU(dimVelocity, SMALL);
        U *= min(scalar(1), fRei*uTau/max(mag(U), maxU));
        correctProcessorPatches(U);
    }

    if (tepsilon.valid())
    {
        tepsilon.ref() = epsilon;
        correctProcessorPatches(tepsilon.ref());
    }

    if (tk.valid())
    {
        tk.ref() = k;
        correctProcessorPatches(tk.ref());
    }

    if (tomega.valid())
    {
        const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);
        tomega.ref() = Cmu*epsilon/(k + k0);
        correctProcessorPatches(tomega.ref());
    }

    if (tR.valid())
    {
        auto& R = tR.ref();

        // (M:Eq. 3)
        const volSphericalTensorField Rdiag(k*twoThirdsI);
        forAll(R, celli)
        {
            R[celli] = Rdiag[celli];
        }
        correctProcessorPatches(R);
    }


    // Write operand fields

    Info<< endl;

    if (initU)
    {
        InfoField(tU->name());
        tU->write();
    }

    if (tepsilon.valid())
    {
        InfoField(tepsilon->name());
        tepsilon->write();
    }

    if (tk.valid())
    {
        InfoField(tk->name());
        tk->write();
    }

    if (tomega.valid())
    {
        InfoField(tomega->name());
        tomega->write();
    }

    if (tR.valid())
    {
        InfoField(tR->name());
        tR->write();
    }

    if (writeF)
    {
        InfoField(f.name());
        f.write();
    }


    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
