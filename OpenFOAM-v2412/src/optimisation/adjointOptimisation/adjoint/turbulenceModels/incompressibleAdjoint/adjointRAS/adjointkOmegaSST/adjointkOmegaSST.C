/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2023 PCOpt/NTUA
    Copyright (C) 2014-2023 FOSS GP
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

#include "adjointkOmegaSST.H"
#include "wallFvPatch.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "linear.H"
#include "reverseLinear.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "omegaWallFunctionFvPatchScalarField.H"
#include "fvOptions.H"
#include "sensitivityTopO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{
namespace adjointRASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointkOmegaSST, 0);
addToRunTimeSelectionTable(adjointRASModel, adjointkOmegaSST, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> adjointkOmegaSST::F1() const
{
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                scalar(500)*nu()/(sqr(y_)*omega())
            ),
            (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> adjointkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
            scalar(500)*nu()/(sqr(y_)*omega())
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


tmp<volScalarField> adjointkOmegaSST::GbyNu
(
    const volScalarField& GbyNu0,
    const volScalarField& F2,
    const volScalarField& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega()
       *max(a1_*omega(), b1_*F2*sqrt(S2))
    );
}


tmp<volScalarField::Internal> adjointkOmegaSST::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega()()
       *max(a1_*omega()(), b1_*F2*sqrt(S2))
    );
}


tmp<volScalarField> adjointkOmegaSST::zeroFirstCell()
{
    auto tzeroFirstCell = volScalarField::New
    (
        "zeroFirstCell",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Foam::one{})
    );
    auto& zeroFirstCell = tzeroFirstCell.ref();

    firstCellIDs_.resize(mesh_.nCells(), -1);
    label counter(0);
    for (const fvPatchScalarField& omegab : omega().boundaryField())
    {
        const fvPatch& patch = omegab.patch();
        if (isA<omegaWallFunctionFvPatchScalarField>(omegab))
        {
            const label patchi = patch.index();
            const labelList& faceCells = patch.faceCells();
            fvPatchScalarField& bf = zeroFirstCell.boundaryFieldRef()[patchi];
            forAll(faceCells, faceI)
            {
                 const label cellI = faceCells[faceI];

                 zeroFirstCell[cellI] = 0.;
                 bf[faceI] = 0.;
                 firstCellIDs_[counter++] = cellI;
            }
        }
    }
    firstCellIDs_.setSize(counter);

    zeroFirstCell.correctBoundaryConditions();

    return tzeroFirstCell;
}


tmp<volScalarField> adjointkOmegaSST::dR_dnut()
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.UaInst();
    word scheme("coeffsDiff");
    tmp<volScalarField> tdR_dnut
    (
        dNutdbMult(U, Ua, nutRef(), scheme)
      + dNutdbMult(k(), ka(), alphaK_, nutRef(), scheme)
      + dNutdbMult(omega(), zeroFirstCell_*wa(), alphaOmega_, nutRef(), scheme)
      - case_1_Pk_*ka()*GbyNu0_*zeroFirstCell_
    );
    volScalarField& dRdnut = tdR_dnut.ref();

    forAll(mesh_.boundary(), pI)
    {
        const fvPatch& patch = mesh_.boundary()[pI];
        const fvPatchScalarField& nutb = nutRef().boundaryField()[pI];
        const labelList& faceCells = patch.faceCells();
        if (isA<nutkWallFunctionFvPatchScalarField>(nutb))
        {
            fvPatchScalarField& bf = dRdnut.boundaryFieldRef()[pI];
            const scalarField kSnGrad(k().boundaryField()[pI].snGrad());
            const scalarField omegaSnGrad
            (
                omega().boundaryField()[pI].snGrad()
            );
            const fvPatchScalarField& akb = alphaK_.boundaryField()[pI];
            const fvPatchScalarField& aOmegab = alphaOmega_.boundaryField()[pI];
            const vectorField USnGrad(U.boundaryField()[pI].snGrad());
            const fvPatchTensorField& gradUb = gradU_.boundaryField()[pI];
            const vectorField nf(mesh_.boundary()[pI].nf());
            forAll(faceCells, fI)
            {
                const label cI(faceCells[fI]);
                bf[fI] =
                  - wa()[cI]*zeroFirstCell_[cI]*aOmegab[fI]*omegaSnGrad[fI]
                  - ka()[cI]*akb[fI]*kSnGrad[fI]
                  - (Ua[cI] & (USnGrad[fI] + (dev2(gradUb[fI]) & nf[fI])));
            }
        }
    }

    return tdR_dnut;
}


tmp<volScalarField> adjointkOmegaSST::dnut_domega
(
    const volScalarField& F2,
    const volScalarField& S,
    const volScalarField& case_1_nut,
    const volScalarField& case_2_nut,
    const volScalarField& case_3_nut
) const
{
    return
        (
          - case_1_nut*k()/sqr(omega())
          - a1_*k()/(b1_*S*F2*F2)*dF2_domega(F2, case_2_nut, case_3_nut)
        );
}


tmp<volScalarField> adjointkOmegaSST::dnut_dk
(
    const volScalarField& F2,
    const volScalarField& S,
    const volScalarField& case_2_nut
) const
{
    return
    (
        a1_/max(a1_*omega(), b1_*F2*S)
      - a1_*k()/(b1_*S*F2*F2)*dF2_dk(F2, case_2_nut)
    );
}


tmp<volScalarField> adjointkOmegaSST::dF2_domega
(
    const volScalarField& F2,
    const volScalarField& case_2_nut,
    const volScalarField& case_3_nut
) const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
            scalar(500)*nu()/(sqr(y_)*omega())
        ),
        scalar(100)
    );

    return
      - scalar(2)*arg2*(1 - F2*F2)*
        (
            case_2_nut*scalar(2)*sqrt(k())/(betaStar_*sqr(omega())*y_)
          + case_3_nut*scalar(500)*nu()/sqr(omega()*y_)
        );
}


tmp<volScalarField> adjointkOmegaSST::dF2_dk
(
    const volScalarField& F2,
    const volScalarField& case_2_nut
) const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
            scalar(500)*nu()/(sqr(y_)*omega())
        ),
        scalar(100)
    );

    return
        case_2_nut*scalar(2)*arg2*(1 - F2*F2)
       /(betaStar_*omega()*y_*sqrt(k()));
}


tmp<volScalarField> adjointkOmegaSST::dGPrime_domega() const
{
    return
        (
            case_2_GPrime_*c1_*betaStar_/a1_
           *(
                max(a1_*omega(), b1_*F2_*S_)
              + case_1_nut_*a1_*omega()
              + (scalar(1) - case_1_nut_)*omega()*b1_*S_
               *dF2_domega(F2_, case_2_nut_, case_3_nut_)
            )
        );
}


tmp<volScalarField> adjointkOmegaSST::dGPrime_dk() const
{
    return
        case_2_GPrime_*c1_*betaStar_/a1_*omega()*b1_*S_
       *dF2_dk(F2_, case_2_nut_);
}


tmp<volScalarField> adjointkOmegaSST::dR_dF1() const
{
    const volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();
    tmp<volScalarField> tGPrime = GbyNu(GbyNu0_, F2_, S2_);
    tmp<volScalarField> tdivU = fvc::div(fvc::absolute(phi, U));
    word scheme("coeffsDiff");

    tmp<volScalarField> tdRdF1
    (
        nutRef()
       *(
            coeffsDifferentiation(k(), ka(), scheme)*(alphaK1_ - alphaK2_)
          + coeffsDifferentiation(omega(), zeroFirstCell_*wa(), scheme)
           *(alphaOmega1_ - alphaOmega2_)
        )
    );
    volScalarField& dRdF1 = tdRdF1.ref();

    typedef fixedValueFvPatchScalarField fixedValue;
    typedef zeroGradientFvPatchScalarField zeroGrad;
    typedef omegaWallFunctionFvPatchScalarField omegaWF;
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& kb = k().boundaryField()[pI];
        const fvPatchScalarField& omegab = omega().boundaryField()[pI];
        fvPatchScalarField& dRdF1b = dRdF1.boundaryFieldRef()[pI];
        if
        (
            isA<zeroGrad>(kb)
         && (isA<zeroGrad>(omegab) || isA<omegaWF>(omegab))
        )
        {
            dRdF1b = dRdF1b.patchInternalField();
        }
        else if (isA<fixedValue>(kb) && isA<fixedValue>(omegab)
        )
        {
            // Note: might need revisiting
            dRdF1b = dRdF1b.patchInternalField();
        }
    }

    volScalarField dR_dF1_coeffs
    (
        zeroFirstCell_*wa()
       *(
            (gamma1_ - gamma2_)*((2.0/3.0)*omega()*tdivU - tGPrime)
          + omega()*omega()*(beta1_ - beta2_) + CDkOmega_
        )
    );

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& kb = k().boundaryField()[pI];
        const fvPatchScalarField& omegab = omega().boundaryField()[pI];
        fvPatchScalarField& dRdF1b = dR_dF1_coeffs.boundaryFieldRef()[pI];
        if
        (
            isA<zeroGrad>(kb)
         && (isA<zeroGrad>(omegab) || isA<omegaWF>(omegab))
        )
        {
            dRdF1b = Zero;
        }
        else if (isA<fixedValue>(kb) && isA<fixedValue>(omegab))
        {
            dRdF1b = dRdF1b.patchInternalField();
        }
    }

    dRdF1 += dR_dF1_coeffs;
    return tdRdF1;
}


tmp<volScalarField> adjointkOmegaSST::dF1_domega
(
    const volScalarField& arg1
) const
{
    return
        (
            4*pow3(arg1)*(scalar(1) - F1_*F1_)
           *(
              - case_1_F1_*sqrt(k())/(betaStar_*omega()*omega()*y_)
              - case_2_F1_*scalar(500)*nu()/sqr(omega()*y_)
              + case_3_F1_*scalar(4)*alphaOmega2_*k()/sqr(CDkOmegaPlus_*y_)
               *CDkOmega_/omega()
            )
        );
}


tmp<volVectorField> adjointkOmegaSST::dF1_dGradOmega
(
    const volScalarField& arg1
) const
{
    return
        (
          - case_3_F1_*scalar(4)*pow3(arg1)*(scalar(1) - F1_*F1_)
           *scalar(8)*k()*sqr(alphaOmega2_/(CDkOmegaPlus_*y_))/omega()*gradK_
        );
}


tmp<volScalarField> adjointkOmegaSST::waEqnSourceFromF1() const
{
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                scalar(500)*nu()/(sqr(y_)*omega())
            ),
            (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
        ),
        scalar(10)
    );

    volScalarField dR_dF1(this->dR_dF1());
    volScalarField dF1_domega(this->dF1_domega(arg1));
    volVectorField dF1_dGradOmega(this->dF1_dGradOmega(arg1));

    surfaceScalarField dR_dGradOmega
    (
        interpolationScheme<vector>("div(dR_dGradOmega)")().
            interpolate(dR_dF1*dF1_dGradOmega) & mesh_.Sf()
    );

    return
        (
            dR_dF1*dF1_domega
          - fvc::div(dR_dGradOmega)
        );
}


tmp<fvScalarMatrix> adjointkOmegaSST::waEqnSourceFromCDkOmega() const
{
    tmp<volVectorField> tsource
    (
        2*zeroFirstCell_*(1 - F1_)*alphaOmega2_/omega()*wa()*gradK_
    );
    volVectorField& source = tsource.ref();

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& omegab = omega().boundaryField()[pI];
        fvPatchVectorField& sourceb = source.boundaryFieldRef()[pI];
        if
        (
            isA<zeroGradientFvPatchScalarField>(omegab)
         || isA<omegaWallFunctionFvPatchScalarField>(omegab)
        )
        {
            sourceb = Zero;
        }
        else if (isA<fixedValueFvPatchScalarField>(omegab))
        {
            sourceb = sourceb.patchInternalField();
        }
    }

    surfaceScalarField sourceFaces
    (
        interpolationScheme<vector>("sourceAdjontkOmegaSST")().
            interpolate(source) & mesh_.Sf()
    );

    return
        (
            // Differentiation of omega in CDkOmega
            fvm::SuSp(zeroFirstCell_*(1. - F1_)*CDkOmega_/omega(), wa())
            // Differentiation of grad(omega) in CDkOmega
          + fvc::div(sourceFaces)
        );
}


tmp<volScalarField> adjointkOmegaSST::kaEqnSourceFromCDkOmega() const
{
    tmp<volVectorField> tsource
    (
        2.*zeroFirstCell_*(1. - F1_)*alphaOmega2_/omega()*wa()*gradOmega_
    );
    volVectorField& source = tsource.ref();

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& kb = k().boundaryField()[pI];
        fvPatchVectorField& sourceb = source.boundaryFieldRef()[pI];
        if (isA<zeroGradientFvPatchScalarField>(kb))
        {
            sourceb = Zero;
        }
        else if (isA<fixedValueFvPatchScalarField>(kb))
        {
            sourceb = sourceb.patchInternalField();
        }
    }

    return
        fvc::div
        (
            interpolationScheme<vector>("sourceAdjontkOmegaSST")().
                interpolate(source) & mesh_.Sf()
        );
}


tmp<volScalarField> adjointkOmegaSST::dF1_dk
(
    const volScalarField& arg1
) const
{
    return
        (
            4*pow3(arg1)*(scalar(1) - F1_*F1_)
           *(
                case_1_F1_*0.5/(betaStar_*omega()*sqrt(k())*y_)
              + case_4_F1_*scalar(4)*alphaOmega2_/(CDkOmegaPlus_*sqr(y_))
            )
        );
}


tmp<volVectorField> adjointkOmegaSST::dF1_dGradK
(
    const volScalarField& arg1
) const
{
    return
        (
          - case_3_F1_*scalar(4)*pow3(arg1)*(scalar(1) - F1_*F1_)
           *scalar(8)*k()*sqr(alphaOmega2_/(CDkOmegaPlus_*y_))/omega()
           *gradOmega_
        );
}


tmp<volScalarField> adjointkOmegaSST::kaEqnSourceFromF1() const
{
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                scalar(500)*nu()/(sqr(y_)*omega())
            ),
            (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
        ),
        scalar(10)
    );

    volScalarField dR_dF1(this->dR_dF1());
    volScalarField dF1_dk(this->dF1_dk(arg1));
    volVectorField dF1_dGradK(this->dF1_dGradK(arg1));

    surfaceScalarField dR_dGradK
    (
        interpolationScheme<vector>("div(dR_dGradK)")().
            interpolate(dR_dF1*dF1_dGradK) & mesh_.Sf()
    );

    return
        (
            dR_dF1*dF1_dk
          - fvc::div(dR_dGradK)
        );
}


tmp<volScalarField> adjointkOmegaSST::coeffsDifferentiation
(
    const volScalarField& primalField,
    const volScalarField& adjointField,
    const word& schemeName
) const
{
    tmp<surfaceInterpolationScheme<scalar>> scheme
        (interpolationScheme<scalar>(schemeName));
    surfaceScalarField snGradPrimal(fvc::snGrad(primalField)*mesh_.magSf());
    surfaceScalarField flux(scheme().interpolate(adjointField)*snGradPrimal);

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& primalbf = primalField.boundaryField()[pI];
        if (isA<zeroGradientFvPatchScalarField>(primalbf))
        {
            // Needless, but just to be safe
            snGradPrimal.boundaryFieldRef()[pI] = Zero;
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchScalarField>(primalbf))
        {
            // Note: to be potentially revisited
            snGradPrimal.boundaryFieldRef()[pI] = Zero;
            flux.boundaryFieldRef()[pI] = Zero;
        }
    }

    return (fvc::div(flux) - adjointField*fvc::div(snGradPrimal));
}


tmp<volScalarField> adjointkOmegaSST::dNutdbMult
(
    const volScalarField& primalField,
    const volScalarField& adjointField,
    const volScalarField& coeffField,
    const volScalarField& bcField,
    const word& schemeName
) const
{
    tmp<surfaceInterpolationScheme<scalar>> scheme
        (interpolationScheme<scalar>(schemeName));
    surfaceScalarField snGradPrimal(fvc::snGrad(primalField)*mesh_.magSf());
    surfaceScalarField flux(scheme().interpolate(adjointField)*snGradPrimal);

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& bc = bcField.boundaryField()[pI];
        if (isA<zeroGradientFvPatchScalarField>(bc))
        {
            const fvPatchScalarField& coeffFieldb =
                coeffField.boundaryField()[pI];
            snGradPrimal.boundaryFieldRef()[pI] *=
                coeffFieldb/coeffFieldb.patchInternalField();
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchScalarField>(bc))
        {
            snGradPrimal.boundaryFieldRef()[pI] = Zero;
            flux.boundaryFieldRef()[pI] = Zero;
        }
    }

    return coeffField*(fvc::div(flux) - adjointField*fvc::div(snGradPrimal));
}


tmp<volScalarField> adjointkOmegaSST::dNutdbMult
(
    const volVectorField& U,
    const volVectorField& Ua,
    const volScalarField& nut,
    const word& schemeName
) const
{
    tmp<surfaceInterpolationScheme<vector>> scheme
        (interpolationScheme<vector>(schemeName));
    surfaceVectorField snGradU(fvc::snGrad(U)*mesh_.magSf());
    surfaceScalarField flux(scheme().interpolate(Ua) & snGradU);

    // Terms form the Laplacian part of the momentum stresses
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& bc = nut.boundaryField()[pI];
        if (isA<zeroGradientFvPatchScalarField>(bc))
        {
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchScalarField>(bc))
        {
            snGradU.boundaryFieldRef()[pI] = Zero;
            flux.boundaryFieldRef()[pI] = Zero;
        }
    }

    // Terms form the tranpose gradient in the momentum stresses
    const surfaceVectorField& Sf = mesh_.Sf();
    surfaceTensorField fluxTranspose
    (
        reverseLinear<vector>(mesh_).interpolate(Ua)*Sf
    );
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& Ub = U.boundaryField()[pI];
        if (!isA<coupledFvPatchVectorField>(Ub))
        {
            const vectorField Uai(Ua.boundaryField()[pI].patchInternalField());
            const vectorField& Sfb = Sf.boundaryField()[pI];
            fluxTranspose.boundaryFieldRef()[pI] = Uai*Sfb;
        }
    }
    volScalarField M(dev2(gradU_) && fvc::div(fluxTranspose));
    const DimensionedField<scalar, volMesh>& V = mesh_.V();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchScalarField& bc = nut.boundaryField()[pI];
        if (isA<zeroGradientFvPatchScalarField>(bc))
        {
            const vectorField Uai(Ua.boundaryField()[pI].patchInternalField());
            const tensorField dev2GradU
            (
                dev2(gradU_.boundaryField()[pI].patchInternalField())
            );
            const vectorField& Sfb = Sf.boundaryField()[pI];
            const labelList& faceCells = mesh_.boundary()[pI].faceCells();
            forAll(faceCells, fI)
            {
                const label celli = faceCells[fI];
                M[celli] -= ((Uai[fI] & dev2GradU[fI]) & Sfb[fI])/V[celli];
            }
        }
    }
    M.correctBoundaryConditions();

  //surfaceScalarField fluxTranspose =
  //    fvc::interpolate(UaGradU, schemeName) & mesh_.Sf();
  //forAll(mesh_.boundary(), pI)
  //{
  //    const fvPatchScalarField& bc = nut.boundaryField()[pI];
  //    const vectorField& Sf = mesh_.boundary()[pI].Sf();
  //    if (isA<zeroGradientFvPatchScalarField>(bc))
  //    {
  //        fluxTranspose.boundaryFieldRef()[pI] =
  //            (
  //                UaGradU.boundaryField()[pI].patchInternalField()
  //              - (
  //                    Ua.boundaryField()[pI].patchInternalField()
  //                  & gradU_.boundaryField()[pI]
  //                )
  //            ) & Sf;
  //    }
  //    else if (isA<fixedValueFvPatchScalarField>(bc))
  //    {
  //        fluxTranspose.boundaryFieldRef()[pI] =
  //            UaGradU.boundaryField()[pI].patchInternalField() & Sf;
  //    }
  //}
    return
        fvc::div(flux) - (Ua & fvc::div(snGradU)) + M;
      //fvc::div(flux + fluxTranspose) - (Ua & fvc::div(snGradU));
}


tmp<volVectorField> adjointkOmegaSST::convectionMeanFlowSource
(
    const volScalarField& primalField,
    const volScalarField& adjointField
) const
{
    // Grab the interpolation scheme from the primal convection term
    tmp<surfaceInterpolationScheme<scalar>> primalInterpolationScheme
    (
        convectionScheme(primalField.name())
    );

    surfaceVectorField interpolatedPrimal
    (
        primalInterpolationScheme().interpolate(primalField)*mesh_.Sf()
    );
    surfaceVectorField flux
    (
      //reverseLinear<scalar>(mesh_).interpolate(adjointField)
        linear<scalar>(mesh_).interpolate(adjointField)
       *interpolatedPrimal
    );

    const volVectorField& U = primalVars_.U();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& bc = U.boundaryField()[pI];
        if (isA<zeroGradientFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchVectorField>(bc))
        {
            interpolatedPrimal.boundaryFieldRef()[pI] = Zero;
            flux.boundaryFieldRef()[pI] = Zero;
        }
    }

    return (-fvc::div(flux) + adjointField*fvc::div(interpolatedPrimal));
}


tmp<volVectorField> adjointkOmegaSST::GMeanFlowSource
(
    tmp<volSymmTensorField>& GbyNuMult
) const
{
    surfaceVectorField flux
    (
        mesh_.Sf() & reverseLinear<symmTensor>(mesh_).interpolate(GbyNuMult())
    );

    const volVectorField& U = primalVars_.U();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& bc = U.boundaryField()[pI];
        if (isA<zeroGradientFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] =
                mesh_.boundary()[pI].Sf()
              & GbyNuMult().boundaryField()[pI].patchInternalField();
        }
    }

    return fvc::div(flux);
}


tmp<volVectorField> adjointkOmegaSST::divUMeanFlowSource
(
    tmp<volScalarField>& divUMult
) const
{
    surfaceVectorField flux
    (
        mesh_.Sf()*reverseLinear<scalar>(mesh_).interpolate(divUMult())
    );

    const volVectorField& U = primalVars_.U();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& bc = U.boundaryField()[pI];
        if (isA<zeroGradientFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] =
                mesh_.boundary()[pI].Sf()
               *divUMult().boundaryField()[pI].patchInternalField();
        }
    }

    divUMult.clear();

    return -fvc::div(flux);
}


tmp<volScalarField> adjointkOmegaSST::diffusionNutMeanFlowMult
(
    const volScalarField& primalField,
    const volScalarField& adjointField,
    const volScalarField& coeffField
) const
{
    // Note: we could grab the snGrad scheme from the diffusion term directly
    surfaceScalarField snGradPrimal(fvc::snGrad(primalField)*mesh_.magSf());
    surfaceScalarField flux
    (
        reverseLinear<scalar>(mesh_).interpolate(adjointField)*snGradPrimal
    );

    const volVectorField& U = primalVars_.U();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& bc = U.boundaryField()[pI];
        if (!isA<coupledFvPatchVectorField>(bc))
        {
            flux.boundaryFieldRef()[pI] = Zero;
            snGradPrimal.boundaryFieldRef()[pI] = Zero;
        }
    }
    return (fvc::div(flux) - adjointField*fvc::div(snGradPrimal))*coeffField;
}


tmp<volVectorField> adjointkOmegaSST::nutMeanFlowSource
(
    tmp<volScalarField>& mult,
    const volScalarField& F2,
    const volScalarField& S,
    const volScalarField& case_1_nut,
    const volTensorField& gradU
) const
{
    volSymmTensorField M
    (
        mult*a1_*k()*(1 - case_1_nut)/(b1_*F2*S*S*S)*twoSymm(gradU)
    );
    M.correctBoundaryConditions();
    mult.clear();

    surfaceVectorField returnFieldFlux
    (
       mesh_.Sf() & reverseLinear<symmTensor>(mesh_).interpolate(M)
    );

    const volVectorField& U = primalVars_.U();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& bc = U.boundaryField()[pI];
        if (isA<zeroGradientFvPatchVectorField>(bc))
        {
            returnFieldFlux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchVectorField>(bc))
        {
            returnFieldFlux.boundaryFieldRef()[pI] =
                mesh_.boundary()[pI].Sf()
              & M.boundaryField()[pI].patchInternalField();
        }
    }

    // Note: potentially missing contributions form patches with a calculated
    // nut bc

    return fvc::div(returnFieldFlux);
}


void adjointkOmegaSST::addWallFunctionTerms
(
    fvScalarMatrix& kaEqn,
    const volScalarField& dR_dnut
)
{
    // Add source term from the differentiation of nutWallFunction
    scalarField& source = kaEqn.source();
    const DimensionedField<scalar, volMesh>& V = mesh_.V();
    const volScalarField& ka = this->ka();
    const volScalarField& wa = this->wa();
    const volScalarField& k = this->k();
    const volScalarField& omega = this->omega();
    forAll(nutRef().boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& nutWall = nutRef().boundaryFieldRef()[patchi];
        if (isA<nutkWallFunctionFvPatchScalarField>(nutWall))
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            const scalarField& magSf = patch.magSf();

            const autoPtr<incompressible::turbulenceModel>& turbModel =
                primalVars_.turbulence();

            const scalarField& y = turbModel().y()[patchi];
            const tmp<scalarField> tnuw = turbModel().nu(patchi);
            const scalarField& nuw = tnuw();

            const nutWallFunctionFvPatchScalarField& nutWF =
                refCast<nutWallFunctionFvPatchScalarField>(nutWall);
            const wallFunctionCoefficients& wallCoeffs = nutWF.wallCoeffs();
            const scalar Cmu = wallCoeffs.Cmu();
            const scalar kappa = wallCoeffs.kappa();
            const scalar E = wallCoeffs.E();
            const scalar yPlusLam = wallCoeffs.yPlusLam();

            const scalar Cmu25 = pow025(Cmu);

            const labelList& faceCells = patch.faceCells();
            const fvPatchScalarField& dR_dnutw =
                dR_dnut.boundaryField()[patchi];
            const fvPatchScalarField& omegaw = omega.boundaryField()[patchi];
            bool addTermsFromOmegaWallFuction
                (isA<omegaWallFunctionFvPatchScalarField>(omegaw));

            const fvPatchVectorField& Uw =
                primalVars_.U().boundaryField()[patchi];
            const scalarField magGradUw(mag(Uw.snGrad()));
            forAll(nuw, facei)
            {
                const label celli = faceCells[facei];

                const scalar sqrtkCell(sqrt(k[celli]));
                const scalar yPlus = Cmu25*y[facei]*sqrtkCell/nuw[facei];
                const scalar logEyPlus = log(E*yPlus);
                const scalar dnut_dyPlus =
                        nuw[facei]*kappa*(logEyPlus - 1)/sqr(logEyPlus);
                const scalar dyPlus_dk =
                    Cmu25*y[facei]/(2*nuw[facei]*sqrtkCell);
                const scalar dnut_dk = dnut_dyPlus*dyPlus_dk;

                if (yPlusLam < yPlus)
                {
                    // Terms from nutkWallFunction
                    source[celli] -= dR_dnutw[facei]*dnut_dk*magSf[facei];
                }
                if (addTermsFromOmegaWallFuction)
                {
                    const scalar denom(Cmu25*kappa*y[facei]);
                    const scalar omegaLog(sqrtkCell/denom);
                    // Terms from omegaLog in omegaWallFunction
                    source[celli] +=
                        wa[celli]*omegaLog/omega[celli]
                       /(2*sqrtkCell*denom);

                    // Terms from G at the first cell off the wall
                    source[celli] +=
                        case_1_Pk_[celli]*ka[celli]*V[celli]
                       *(
                            (nutWall[facei] + nuw[facei])
                           *magGradUw[facei]
                           *Cmu25
                           /(2.0*sqrtkCell*kappa*y[facei])
                        );

                    if (yPlusLam < yPlus)
                    {
                        source[celli] +=
                            case_1_Pk_[celli]*ka[celli]*V[celli]
                           *dnut_dk
                           *magGradUw[facei]
                           *Cmu25*sqrtkCell
                           /(kappa*y[facei]);
                    }
                }
            }
        }
    }
}


void adjointkOmegaSST::updatePrimalRelatedFields()
{
    if (changedPrimalSolution_)
    {
        Info<< "Updating primal-based fields of the adjoint turbulence "
            << "model ..." << endl;

        const volVectorField& U = primalVars_.U();
        gradU_ = fvc::grad(U);
        gradOmega_ = fvc::grad(omega());
        gradK_ = fvc::grad(k());

        S2_ = 2*magSqr(symm(gradU_))
            + dimensionedScalar(dimless/sqr(dimTime), 1.e-21);
        S_ = sqrt(S2_);
        GbyNu0_ = gradU_ && devTwoSymm(gradU_);

        // Instead of computing G directly here, delegate to RASModelVariables
        // to make sure G is computed consistently when the primal fields are
        // averaged. The local value computed here is just used to update
        // the switch fields
        /*
        volScalarField G(primalVars_.turbulence()->GName(), nutRef()*GbyNu0_);
        omega().correctBoundaryConditions();
        */
        volScalarField G
        (
            IOobject
            (
                IOobject::scopedName(type(), "G"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimArea/pow3(dimTime), Zero)
        );
        G.ref() = primalVars_.RASModelVariables()->G()();

        CDkOmega_ =
            (2*alphaOmega2_)*(gradK_ & gradOmega_)/omega();
        CDkOmegaPlus_ = max
        (
            CDkOmega_,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        );
        F1_ = F1();
        F2_ = F2();
        case_1_Pk_ = neg(G - c1_*betaStar_*k()*omega());
        case_2_Pk_ = pos0(G - c1_*betaStar_*k()*omega());
        case_3_Pk_ = neg(a1_*omega() - b1_*F2_*S_);


        alphaK_ = alphaK(F1_);
        alphaOmega_ = alphaOmega(F1_);
        beta_ = beta(F1_);
        gamma_ = gamma(F1_);

        // Switch fields for F1
        {
            volScalarField arg1_C3
            (
                (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_)
              - scalar(500)*nu()/(sqr(y_)*omega())
            );
            volScalarField arg1_C2
            (
                max
                (
                    (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                    scalar(500)*nu()/(sqr(y_)*omega())
                )
              - (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
            );
            volScalarField arg1_C1
            (
                min
                (
                    max
                    (
                        (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                        scalar(500)*nu()/(sqr(y_)*omega())
                    ),
                    (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
                ) - scalar(10)
            );
            volScalarField CDkOmegaPlus_C1
            (
                CDkOmega_
              - dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
            );

            case_1_F1_ = pos(arg1_C3)*neg(arg1_C2)*neg(arg1_C1);
            case_2_F1_ = neg0(arg1_C3)*neg(arg1_C2)*neg(arg1_C1);
            case_3_F1_ = pos0(arg1_C2)*neg(arg1_C1)*pos(CDkOmegaPlus_C1);
            case_4_F1_ = pos0(arg1_C2)*neg(arg1_C1);
        }

        // Switch fields for nut
        {
            volScalarField nut_C1(a1_*omega() - b1_*F2_*S_);
            volScalarField arg2_C2
            (
                (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_)
              - scalar(500)*nu()/(sqr(y_)*omega())
            );
            volScalarField arg2_C1
            (
                max
                (
                    (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
                    scalar(500)*nu()/(sqr(y_)*omega())
                ) - scalar(100)
            );

            case_1_nut_ = pos(nut_C1);
            case_2_nut_ = neg0(nut_C1)*pos(arg2_C2)*neg(arg2_C1);
            case_3_nut_ = neg0(nut_C1)*neg0(arg2_C2)*neg(arg2_C1);
        }

        {
            volScalarField GPrime_C1
            (
                GbyNu0_
              - (c1_/a1_)*betaStar_*omega()*max(a1_*omega(), b1_*F2_*S_)
            );
            case_1_GPrime_ = neg(GPrime_C1);
            case_2_GPrime_ = pos0(GPrime_C1);
        }

        dnut_domega_ =
            dnut_domega(F2_, S_, case_1_nut_, case_2_nut_, case_3_nut_);
        dnut_dk_ = dnut_dk(F2_, S_, case_2_nut_);
        DOmegaEff_ = DomegaEff(F1_);
        DkEff_ = DkEff(F1_);

        changedPrimalSolution_ = false;
    }
}


tmp<surfaceInterpolationScheme<scalar>> adjointkOmegaSST::convectionScheme
(
    const word& varName
) const
{
    const surfaceScalarField& phi = primalVars_.phi();
    const surfaceScalarField& phiInst = primalVars_.phiInst();
    word divEntry("div(" + phiInst.name() + ',' + varName +')');
    ITstream& divScheme = mesh_.divScheme(divEntry);
    // Skip the first entry which might be 'bounded' or 'Gauss'.
    // If it is 'bounded', skip the second entry as well
    word discarded(divScheme);
    if (discarded == "bounded")
    {
        discarded = word(divScheme);
    }
    return surfaceInterpolationScheme<scalar>::New(mesh_, phi, divScheme);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointkOmegaSST::adjointkOmegaSST
(
    incompressibleVars& primalVars,
    incompressibleAdjointMeanFlowVars& adjointVars,
    objectiveManager& objManager,
    const word& adjointTurbulenceModelName,
    const word& modelName
)
:
    adjointRASModel
    (
        modelName,
        primalVars,
        adjointVars,
        objManager,
        adjointTurbulenceModelName
    ),

    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(primalVars_.RASModelVariables()().d()),

    //Primal Gradient Fields
    gradU_
    (
        IOobject
        (
            "rasModel::gradU",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimless/dimTime, Zero)
    ),
    gradOmega_
    (
        IOobject
        (
            "rasModel::gradOmega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless/dimTime/dimLength, Zero)
    ),
    gradK_
    (
        IOobject
        (
            "rasModel::gradK",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/sqr(dimTime), Zero)
    ),

    S2_
    (
        IOobject
        (
            "S2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/sqr(dimTime), Zero)
    ),
    S_
    (
        IOobject
        (
            "kOmegaSST_S",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),
    GbyNu0_
    (
        IOobject
        (
            "adjointRASModel::GbyNu0",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/sqr(dimTime), Zero)
    ),
    CDkOmega_
    (
        IOobject
        (
            "CDkOmega_",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/sqr(dimTime), Zero)
    ),
    CDkOmegaPlus_
    (
        IOobject
        (
            "CDkOmegaPlus",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/sqr(dimTime), Zero)
    ),
    F1_
    (
        IOobject
        (
            "F1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    F2_
    (
        IOobject
        (
            "F2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    // Model Field coefficients
    alphaK_
    (
        IOobject
        (
            "alphaK",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    alphaOmega_
    (
        IOobject
        (
            "alphaOmega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    beta_
    (
        IOobject
        (
            "beta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    case_1_F1_
    (
        IOobject
        (
            "case_1_F1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_2_F1_
    (
        IOobject
        (
            "case_2_F1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_3_F1_
    (
        IOobject
        (
            "case_3_F1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_4_F1_
    (
        IOobject
        (
            "case_4_F1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_1_Pk_
    (
        IOobject
        (
            "case_1_Pk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_2_Pk_
    (
        IOobject
        (
            "case_2_Pk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_3_Pk_
    (
        IOobject
        (
            "case_3_Pk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    case_1_nut_
    (
        IOobject
        (
            "case_1_nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_2_nut_
    (
        IOobject
        (
            "case_2_nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_3_nut_
    (
        IOobject
        (
            "case_3_nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_1_GPrime_
    (
        IOobject
        (
            "case_1_GPrime",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    case_2_GPrime_
    (
        IOobject
        (
            "case_2_GPrime",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    // Zero 1rst cell field
    firstCellIDs_(0),
    zeroFirstCell_(zeroFirstCell()),

    // Turbulence model multipliers
    dnut_domega_
    (
        IOobject
        (
            "dnut_domega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength), Zero)
    ),
    dnut_dk_
    (
        IOobject
        (
            "dnut_dk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, Zero)
    ),
    DOmegaEff_
    (
        IOobject
        (
            "DomegaEff",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(nutRef().dimensions(), Zero)
    ),
    DkEff_
    (
        IOobject
        (
            "DkEff",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(nutRef().dimensions(), Zero)
    )
{
    adjointTMVariablesBaseNames_.setSize(2);
    adjointTMVariablesBaseNames_[0] = "ka";
    adjointTMVariablesBaseNames_[1] = "wa";
    // Read in adjoint fields
    variablesSet::setField
    (
        adjointTMVariable1Ptr_,
        mesh_,
        "ka",
        adjointVars.solverName(),
        adjointVars.useSolverNameForFields()
    );
    variablesSet::setField
    (
        adjointTMVariable2Ptr_,
        mesh_,
        "wa",
        adjointVars.solverName(),
        adjointVars.useSolverNameForFields()
    );

    setMeanFields();

    // No sensitivity contributions from the adjoint to the eikonal equation
    // for the moment
    includeDistance_ = false;

    // Update the primal related fields here so that functions computing
    // sensitivities have the updated fields in case of continuation
    updatePrimalRelatedFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> adjointkOmegaSST::devReff() const
{
    const volVectorField& Ua = adjointVars_.UaInst();
    return devReff(Ua);
}


tmp<volSymmTensorField> adjointkOmegaSST::devReff
(
    const volVectorField& U
) const
{
    return tmp<volSymmTensorField>::New
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*devTwoSymm(fvc::grad(U))
        );
}


tmp<fvVectorMatrix> adjointkOmegaSST::divDevReff(volVectorField& Ua) const
{
    tmp<volScalarField> tnuEff = nuEff();
    const volScalarField& nuEff = tnuEff();

    return
    (
      - fvm::laplacian(nuEff, Ua)
      - fvc::div(nuEff*dev(fvc::grad(Ua)().T()))
    );

    /* WIP
    const volVectorField& U = primalVars_.U();
    const surfaceVectorField& Sf = mesh_.Sf();
    tmp<surfaceTensorField> tflux =
        reverseLinear<vector>(mesh_).interpolate(Ua)*Sf;
    surfaceTensorField& flux = tflux.ref();
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& Ub = U.boundaryField()[pI];
        if (!isA<coupledFvPatchVectorField>(Ub))
        {
            const vectorField Uai = Ua.boundaryField()[pI].patchInternalField();
            const vectorField& Sfb = Sf.boundaryField()[pI];
            flux.boundaryFieldRef()[pI] = Uai*Sfb;
        }
    }
    volTensorField M(nuEff*dev2(fvc::div(flux)));
    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& Ub = U.boundaryField()[pI];
        if (!isA<coupledFvPatchVectorField>(Ub))
        {
            const fvPatchScalarField& nuEffb = nuEff.boundaryField()[pI];
            const vectorField nf = mesh_.boundary()[pI].nf();
            const vectorField Uai = Ua.boundaryField()[pI].patchInternalField();
            const labelList& faceCells = mesh_.boundary()[pI].faceCells();
            const vectorField& Sfb = Sf.boundaryField()[pI];

            forAll(faceCells, fI)
            {
                const label celli = faceCells[fI];
                const tensor t(dev2(Uai[fI]*Sfb[fI]));
                M[celli] -= nuEffb[fI]*(t - nf[fI]*(nf[fI] & t))/V[celli];
            }
        }
    }
    M.correctBoundaryConditions();

    surfaceVectorField returnFlux =
      - (Sf & reverseLinear<tensor>(mesh_).interpolate(M));
    forAll(mesh_.boundary(), pI)
    {
        const fvPatchVectorField& Ub = U.boundaryField()[pI];
        if (isA<zeroGradientFvPatchVectorField>(Ub))
        {
            returnFlux.boundaryFieldRef()[pI] = Zero;
        }
        else if (isA<fixedValueFvPatchVectorField>(Ub))
        {
            const scalarField& deltaCoeffs = mesh_.boundary()[pI].deltaCoeffs();
            const fvPatchScalarField& nuEffb = nuEff.boundaryField()[pI];
            const vectorField Uai = Ua.boundaryField()[pI].patchInternalField();
            const vectorField nf = mesh_.boundary()[pI].nf();
            const vectorField& Sfb = Sf.boundaryField()[pI];

            returnFlux.boundaryFieldRef()[pI] =
              - (Sfb & M.boundaryField()[pI].patchInternalField())
              + nuEffb*deltaCoeffs*(nf & dev2(Uai*Sfb));
        }
    }

    return
    (
      - fvm::laplacian(nuEff, Ua)
      + fvc::div(returnFlux)
    );
    */
}


tmp<volVectorField> adjointkOmegaSST::nonConservativeMomentumSource() const
{
    return (ka()*gradK_ + wa()*gradOmega_);
}


tmp<volVectorField> adjointkOmegaSST::adjointMeanFlowSource()
{
    tmp<volVectorField> tmeanFlowSource
    (
        tmp<volVectorField>::New
        (
            IOobject
            (
                "adjointMeanFlowSource" + type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimVelocity/dimTime, Zero)
        )
    );
    volVectorField& meanFlowSource = tmeanFlowSource.ref();

    // Contributions from the convection terms of the turbulence model
    meanFlowSource +=
        convectionMeanFlowSource(omega(), zeroFirstCell_*wa())
      + convectionMeanFlowSource(k(), ka());

    // Contributions from GbyNu, including gradU
    tmp<volSymmTensorField> twoSymmGradU(twoSymm(gradU_));
    tmp<volSymmTensorField> GbyNuMult
    (
        // First part of GPrime and G from Pk
        2.*dev(twoSymmGradU())*zeroFirstCell_
       *(wa()*gamma_*case_1_GPrime_ + ka()*nutRef()*case_1_Pk_)
        // Second part of GPrime
      + twoSymmGradU()*wa()*zeroFirstCell_*gamma_*case_2_GPrime_
       *(1. - case_1_nut_)*c1_/a1_*betaStar_*omega()*b1_*F2_/S_
    );
    twoSymmGradU.clear();
    meanFlowSource += GMeanFlowSource(GbyNuMult);
    GbyNuMult.clear();

    // Contributions from divU
    tmp<volScalarField> divUMult
    (
        (2.0/3.0)*(zeroFirstCell_*wa()*omega()*gamma_ + ka()*k())
    );
    meanFlowSource += divUMeanFlowSource(divUMult);

    // Contributions from S2, existing in nut
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.UaInst();
    tmp<volScalarField> nutMeanFlowSourceMult
    (
        // nut in the diffusion coefficients
        diffusionNutMeanFlowMult(k(), ka(), alphaK_)
      + diffusionNutMeanFlowMult(omega(), zeroFirstCell_*wa(), alphaOmega_)
      + dNutdbMult(U, Ua, nutRef(), "coeffsDiff")
        // nut in G
      - ka()*case_1_Pk_*zeroFirstCell_*GbyNu0_
    );
    meanFlowSource +=
        nutMeanFlowSource
        (
            nutMeanFlowSourceMult,
            F2_,
            S_,
            case_1_nut_,
            gradU_
        );

    // G at the first cell includes mag(U.snGrad())
    // Add term here
    forAll(omega().boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& omegaWall = omega().boundaryFieldRef()[patchi];
        if (isA<omegaWallFunctionFvPatchScalarField>(omegaWall))
        {
            const fvPatch& patch = mesh_.boundary()[patchi];

            const autoPtr<incompressible::turbulenceModel>& turbModel =
                primalVars_.turbulence();

            const scalarField& y = turbModel().y()[patchi];
            const tmp<scalarField> tnuw = turbModel().nu(patchi);
            const scalarField& nuw = tnuw();

            const nutWallFunctionFvPatchScalarField& nutw =
                refCast<nutWallFunctionFvPatchScalarField>
                    (nutRef().boundaryFieldRef()[patchi]);
            const wallFunctionCoefficients& wallCoeffs = nutw.wallCoeffs();
            const scalar Cmu = wallCoeffs.Cmu();
            const scalar kappa = wallCoeffs.kappa();
            const scalar Cmu25 = pow025(Cmu);

            const labelList& faceCells = patch.faceCells();

            const fvPatchVectorField& Uw =
                primalVars_.U().boundaryField()[patchi];
            vectorField snGradUw(Uw.snGrad());
            const scalarField& deltaCoeffs = patch.deltaCoeffs();
            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                // Volume will be added when meanFlowSource is added to UaEqn
                meanFlowSource[celli] +=
                    ka()[celli]*case_1_Pk_[celli]
                   *(nutw[facei] + nuw[facei])
                   *snGradUw[facei].normalise()
                   *Cmu25*sqrt(k()[celli])
                   *deltaCoeffs[facei]
                   /(kappa*y[facei]);
            }
        }
    }
    return tmeanFlowSource;
}


tmp<volScalarField> adjointkOmegaSST::nutJacobianTMVar1() const
{
    // Compute dnut_dk anew since the current copy of dnut_dk_
    // might not be updated
    const volVectorField& U = primalVars_.U();
    tmp<volScalarField> S2
    (
        2*magSqr(symm(fvc::grad(U)))
      + dimensionedScalar(dimless/sqr(dimTime), 1.e-21)
    );
    volScalarField S(sqrt(S2));
    volScalarField F2(this->F2());

    // Computation of nut switches
    volScalarField nut_C1(a1_*omega() - b1_*F2*S);
    volScalarField arg2_C2
    (
        (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_)
      - scalar(500)*nu()/(sqr(y_)*omega())
    );
    volScalarField arg2_C1
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
            scalar(500)*nu()/(sqr(y_)*omega())
        ) - scalar(100)
    );
    volScalarField case_2_nut(neg0(nut_C1)*pos(arg2_C2)*neg(arg2_C1));

    return dnut_dk(F2, S, case_2_nut);
}


tmp<volScalarField> adjointkOmegaSST::nutJacobianTMVar2() const
{
    // Compute dnut_omega anew since the current copy of dnut_domega_
    // might not be updated
    const volVectorField& U = primalVars_.U();
    tmp<volScalarField> S2
    (
        2*magSqr(symm(fvc::grad(U)))
      + dimensionedScalar(dimless/sqr(dimTime), 1.e-21)
    );
    volScalarField S(sqrt(S2));
    volScalarField F2(this->F2());

    // Computation of nut switches
    volScalarField nut_C1(a1_*omega() - b1_*F2*S);
    volScalarField arg2_C2
    (
        (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_)
      - scalar(500)*nu()/(sqr(y_)*omega())
    );
    volScalarField arg2_C1
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k())/(omega()*y_),
            scalar(500)*nu()/(sqr(y_)*omega())
        ) - scalar(100)
    );
    volScalarField case_1_nut(pos(nut_C1));
    volScalarField case_2_nut(neg0(nut_C1)*pos(arg2_C2)*neg(arg2_C1));
    volScalarField case_3_nut(neg0(nut_C1)*neg0(arg2_C2)*neg(arg2_C1));

    return dnut_domega(F2, S, case_1_nut, case_2_nut, case_3_nut);
}


tmp<volVectorField> adjointkOmegaSST::nutJacobianU
(
    tmp<volScalarField>& dNutdUMult
) const
{
    const volVectorField& U = primalVars_.U();
    volTensorField gradU(fvc::grad(U));
    tmp<volScalarField> S2
    (
        2*magSqr(symm(gradU))
      + dimensionedScalar(dimless/sqr(dimTime), 1.e-21)
    );
    volScalarField S(sqrt(S2));
    volScalarField F2(this->F2());

    // Computation of nut switches
    volScalarField nut_C1(a1_*omega() - b1_*F2*S);
    volScalarField case_1_nut(pos(nut_C1));

    return nutMeanFlowSource(dNutdUMult, F2, S, case_1_nut, gradU);
}


tmp<scalarField> adjointkOmegaSST::diffusionCoeffVar1(label patchI) const
{
    return
        (
            alphaK_.boundaryField()[patchI]*nutRef().boundaryField()[patchI]
          + nu()().boundaryField()[patchI]
        );
}


tmp<scalarField> adjointkOmegaSST::diffusionCoeffVar2(label patchI) const
{
    return
        (
            alphaOmega_.boundaryField()[patchI]*nutRef().boundaryField()[patchI]
          + nu()().boundaryField()[patchI]
        );
}


void adjointkOmegaSST::correct()
{
    adjointRASModel::correct();

    if (!adjointTurbulence_)
    {
        return;
    }

    updatePrimalRelatedFields();

    // Primal and adjoint fields
    const volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();
    volScalarField dR_dnut(this->dR_dnut());
    volScalarField::Internal divU(fvc::div(fvc::absolute(phi, U)));

    fv::options& fvOptions(fv::options::New(mesh_));

    tmp<fvScalarMatrix> waEqn
    (
        fvm::div(-phi, wa())
      + fvm::SuSp(zeroFirstCell_*fvc::div(phi), wa())
      - fvm::laplacian(DOmegaEff_, wa())
      + waEqnSourceFromCDkOmega()
      + waEqnSourceFromF1()
      + dR_dnut*dnut_domega_
      + fvm::Sp(zeroFirstCell_*scalar(2.)*beta_*omega(), wa())
      + fvm::SuSp
        (
            zeroFirstCell_()*gamma_*((2.0/3.0)*divU - dGPrime_domega().ref()()),
            wa()
        )
      - (case_2_Pk_*c1_ - scalar(1))*betaStar_*k()*ka()
     ==
        fvOptions(wa())
    );

    // Boundary manipulate changes the diagonal component, so relax has to
    // come after that
    waEqn.ref().boundaryManipulate(wa().boundaryFieldRef());
    waEqn.ref().relax();

    // Sources from the objective should be added after the boundary
    // manipulation
    objectiveManager_.addSource(waEqn.ref());
    fvOptions.constrain(waEqn.ref());
    waEqn.ref().solve();

    // Adjoint Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kaEqn
    (
        fvm::div(-phi, ka())
      + fvm::SuSp(fvc::div(phi), ka())
      - fvm::laplacian(DkEff_, ka())
      + fvm::Sp(betaStar_*omega(), ka())
      - case_2_Pk_()*c1_*betaStar_*omega()()*ka()()
      + fvm::SuSp(scalar(2.0/3.0)*divU, ka())
      + kaEqnSourceFromCDkOmega()
      + kaEqnSourceFromF1()
      + dR_dnut()*dnut_dk_()
      - zeroFirstCell_()*gamma_*dGPrime_dk().ref()()*wa()
     ==
        fvOptions(ka())
    );

    kaEqn.ref().relax();
    kaEqn.ref().boundaryManipulate(ka().boundaryFieldRef());
    addWallFunctionTerms(kaEqn.ref(), dR_dnut);
    fvOptions.constrain(kaEqn.ref());
    // Add sources from the objective functions
    objectiveManager_.addSource(kaEqn.ref());

    kaEqn.ref().solve();

    if (adjointVars_.getSolverControl().printMaxMags())
    {
        dimensionedScalar maxwa = max(mag(wa()));
        dimensionedScalar maxka = max(mag(ka()));
        Info<< "Max mag (" << wa().name() << ") = " << maxwa.value() << endl;
        Info<< "Max mag (" << ka().name() << ") = " << maxka.value() << endl;
    }
}


const boundaryVectorField& adjointkOmegaSST::adjointMomentumBCSource() const
{
    return adjMomentumBCSourcePtr_();
}


const boundaryVectorField& adjointkOmegaSST::wallShapeSensitivities()
{
    boundaryVectorField& wallShapeSens = wallShapeSensitivitiesPtr_();
    volTensorField FITerm(FISensitivityTerm());

    forAll(mesh_.boundary(), patchi)
    {
       vectorField nf(mesh_.boundary()[patchi].nf());
       wallShapeSens[patchi] = nf & FITerm.boundaryField()[patchi];
    }
    return wallShapeSens;
}


const boundaryVectorField& adjointkOmegaSST::wallFloCoSensitivities()
{
    return wallFloCoSensitivitiesPtr_();
}


tmp<volScalarField> adjointkOmegaSST::distanceSensitivities()
{
    return tmp<volScalarField>::New
        (
            IOobject
            (
                "adjointEikonalSource" + type(),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength/pow3(dimTime), Zero)
       );
}


tmp<volTensorField> adjointkOmegaSST::FISensitivityTerm()
{
    const volVectorField& U = primalVars_.U();
    const volScalarField& kInst =
        primalVars_.RASModelVariables()->TMVar1Inst();
    const volScalarField& omegaInst =
        primalVars_.RASModelVariables()->TMVar2Inst();

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k())/(omega()*y_),
                scalar(500)*nu()/(sqr(y_)*omega())
            ),
            (4*alphaOmega2_)*k()/(CDkOmegaPlus_*sqr(y_))
        ),
        scalar(10)
    );

    // Interpolation schemes used by the primal convection terms
    auto kScheme(convectionScheme(kInst.name()));
    auto omegaScheme(convectionScheme(omegaInst.name()));
    const surfaceVectorField& Sf = mesh_.Sf();
    tmp<volTensorField> tFISens
    (
        tmp<volTensorField>::New
        (
            IOobject
            (
                type() + "FISensTerm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    volTensorField& FISens = tFISens.ref();
    FISens =
        // k convection
      - ka()*fvc::div
        (
            kScheme().interpolate(k())
           *linear<vector>(mesh_).interpolate(U)*Sf
        )
        // k diffusion
      + ka()*T(fvc::grad(DkEff_*gradK_))
      - DkEff_*(fvc::grad(ka())*gradK_)
        // omega convection
      - wa()*zeroFirstCell_*fvc::div
        (
            omegaScheme().interpolate(omega())
           *linear<vector>(mesh_).interpolate(U)*Sf
        )
        // omega diffusion
      + wa()*zeroFirstCell_*T(fvc::grad(DOmegaEff_*gradOmega_))
      - DOmegaEff_*(fvc::grad(wa()*zeroFirstCell_)*gradOmega_)
        // terms including GbyNu0
      + (
            case_1_GPrime_*wa()*gamma_
          + case_1_Pk_*ka()*nutRef()
        )*2.*T(gradU_ & devTwoSymm(gradU_))*zeroFirstCell_
        // S2 (includes contribution from nut in UEqn as well)
      + (
            dR_dnut()*a1_*k()/(b1_*S_*S_*S_*F2_)
          + wa()*zeroFirstCell_*gamma_*case_2_GPrime_
            *(c1_/a1_)*betaStar_*omega()*b1_*F2_/S_
        )*T(gradU_ & twoSymm(gradU_))*(1. - case_1_nut_)
        // CDkOmega in omegaEqn
      + 2.*wa()*(1. - F1_)*alphaOmega2_/omega()*zeroFirstCell_
        *(gradOmega_*gradK_ + gradK_*gradOmega_)
        // F1
      - dR_dF1()
        *(dF1_dGradK(arg1)*gradK_ + dF1_dGradOmega(arg1)*gradOmega_);

    FISens.correctBoundaryConditions();

    return tFISens;
}


tmp<scalarField> adjointkOmegaSST::topologySensitivities
(
    const word& designVarsName
) const
{
    fv::options& fvOptions(fv::options::New(this->mesh_));
    auto tres(tmp<scalarField>::New(mesh_.nCells(), Zero));

    // Sensitivity from the source term in the k equation
    scalarField auxSens
        (k().primitiveField()*ka().primitiveField());
    sensitivityTopO::postProcessSens
    (
        tres.ref(), auxSens, fvOptions, k().name(), designVarsName
    );

    // Sensitivity from the source term in the omega equation
    auxSens = omega().primitiveField()*wa().primitiveField();
    sensitivityTopO::postProcessSens
    (
        tres.ref(), auxSens, fvOptions, omega().name(), designVarsName
    );

    return tres;
}


void adjointkOmegaSST::nullify()
{
    variablesSet::nullifyField(ka());
    variablesSet::nullifyField(wa());
}


bool adjointkOmegaSST::read()
{
    if (adjointRASModel::read())
    {
        kappa_.readIfPresent(coeffDict());
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adjointRASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
