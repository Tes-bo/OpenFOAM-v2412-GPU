/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "fv.H"
#include "HashTable.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "cyclicAMIFvPatch.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<ddtScheme<Type>> ddtScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing ddtScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Ddt scheme not specified" << endl << endl
            << "Valid ddt schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto* ctorPtr = IstreamConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "ddt",
            schemeName,
            *IstreamConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> ddtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    NotImplemented;
    return nullptr;
}


template<class Type>
tmp<fvMatrix<Type>> ddtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    NotImplemented;
    return nullptr;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> ddtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    NotImplemented;
    return nullptr;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr
)
{
    if (fv::debug)
    {
        InfoInFunction << "Using standard version" << endl;
    }

    tmp<surfaceScalarField> tddtCouplingCoeff
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddtCouplingCoeff",
                U.mesh().time().timeName(),
                U.mesh().thisDb()
            ),
            U.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff.ref();

    if (ddtPhiCoeff_ < 0)
    {
        // v1712 and earlier
        ddtCouplingCoeff -= min
        (
            mag(phiCorr)
           /(mag(phi) + dimensionedScalar("small", phi.dimensions(), SMALL)),
            scalar(1)
        );
    }
    else
    {
        ddtCouplingCoeff =
            dimensionedScalar("ddtPhiCoeff", dimless, ddtPhiCoeff_);
    }

    surfaceScalarField::Boundary& ccbf = ddtCouplingCoeff.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
            U.boundaryField()[patchi].fixesValue()
         || isA<cyclicAMIFvPatch>(mesh().boundary()[patchi])
        )
        {
            ccbf[patchi] = 0.0;
        }
    }

    if (debug > 1)
    {
        InfoInFunction
            << "ddtCouplingCoeff mean max min = "
            << gAverage(ddtCouplingCoeff.primitiveField())
            << " " << gMax(ddtCouplingCoeff.primitiveField())
            << " " << gMin(ddtCouplingCoeff.primitiveField())
            << endl;
    }

    return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeffExperimental
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr
)
{
    if (fv::debug)
    {
        InfoInFunction << "Using experimental version" << endl;
    }

    tmp<surfaceScalarField> tddtCouplingCoeff
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ddtCouplingCoeff",
                U.mesh().time().timeName(),
                U.mesh().thisDb()
            ),
            U.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff.ref();

    if (ddtPhiCoeff_ < 0)
    {
        // See note below re: commented code
        ddtCouplingCoeff -= min
        (
            //  mag(phiCorr)
            // *mesh().time().deltaT()*mag(mesh().deltaCoeffs())/mesh().magSf(),
            //  scalar(1)
            mag(phiCorr)
            *mesh().time().deltaT()*mesh().deltaCoeffs()/mesh().magSf(),
            scalar(1)
        );

        // Note: setting oriented to false to avoid having to use mag(deltaCoeffs)
        // - the deltaCoeffs field is always positive (scalars)
        ddtCouplingCoeff.setOriented(false);
    }
    else
    {
        ddtCouplingCoeff =
            dimensionedScalar("ddtPhiCoeff", dimless, ddtPhiCoeff_);
    }

    surfaceScalarField::Boundary& ccbf = ddtCouplingCoeff.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
            U.boundaryField()[patchi].fixesValue()
         || isA<cyclicAMIFvPatch>(mesh().boundary()[patchi])
        )
        {
            ccbf[patchi] = 0.0;
        }
    }

    if (debug > 1)
    {
        InfoInFunction
            << "ddtCouplingCoeff mean max min = "
            << gAverage(ddtCouplingCoeff.primitiveField())
            << " " << gMax(ddtCouplingCoeff.primitiveField())
            << " " << gMin(ddtCouplingCoeff.primitiveField())
            << endl;
    }

    return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr,
    const volScalarField& rho
)
{
    if (experimentalDdtCorr)
    {
        return
            fvcDdtPhiCoeffExperimental
            (
                U,
                phi,
                phiCorr/fvc::interpolate(rho)
            );
    }
    else
    {
        return fvcDdtPhiCoeff(U, phi, phiCorr);
    }
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    if (experimentalDdtCorr)
    {
        return
            fvcDdtPhiCoeffExperimental
            (
                U,
                phi,
                phi - fvc::dotInterpolate(mesh().Sf(), U)
            );
    }
    else
    {
        return
            fvcDdtPhiCoeff
            (
                U,
                phi,
                phi - fvc::dotInterpolate(mesh().Sf(), U)
            );
    }
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& rhoU,
    const fluxFieldType& phi,
    const volScalarField& rho
)
{
    if (experimentalDdtCorr)
    {
        return fvcDdtPhiCoeffExperimental
        (
            rhoU,
            phi,
            (phi - fvc::dotInterpolate(mesh().Sf(), rhoU))
           /fvc::interpolate(rho)
        );
    }
    else
    {
        return fvcDdtPhiCoeff
        (
            rhoU,
            phi,
            (phi - fvc::dotInterpolate(mesh().Sf(), rhoU))
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
