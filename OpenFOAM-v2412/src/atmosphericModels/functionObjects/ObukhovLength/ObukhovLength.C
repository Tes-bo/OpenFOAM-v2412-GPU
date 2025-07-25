/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020,2023 OpenCFD Ltd
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

#include "ObukhovLength.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "coupledFvPatchField.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ObukhovLength, 0);
    addToRunTimeSelectionTable(functionObject, ObukhovLength, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::functionObjects::ObukhovLength::calcOL()
{
    const auto* rhoPtr = mesh_.findObject<volScalarField>("rho");
    const volScalarField& nut = mesh_.lookupObject<volScalarField>("nut");
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& alphat = mesh_.lookupObject<volScalarField>("alphat");
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    volScalarField* result1 = mesh_.getObjectPtr<volScalarField>(resultName1_);
    volScalarField* result2 = mesh_.getObjectPtr<volScalarField>(resultName2_);

    const bool isNew = !result1;

    if (!result1)
    {
        result1 = new volScalarField
        (
            IOobject
            (
                resultName1_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimLength,
            fvPatchFieldBase::zeroGradientType()
        );

        result1->store();
    }

    if (!result2)
    {
        result2 = new volScalarField
        (
            IOobject
            (
                resultName2_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimVelocity,
            fvPatchFieldBase::zeroGradientType()
        );

        result2->store();
    }

    volScalarField B(alphat*beta_*(fvc::grad(T) & g_));
    if (rhoPtr)
    {
        const auto& rho = *rhoPtr;
        B /= rho;
    }
    else
    {
        const dimensionedScalar rho(dimless, rhoRef_);
        B /= rho;
    }

    *result2 = // Ustar
        sqrt
        (
            max
            (
                nut*sqrt(2*magSqr(symm(fvc::grad(U)))),
                dimensionedScalar(sqr(U.dimensions()), VSMALL)
            )
        );


    // Special bit of coding here to handle cyclics
    // - sign(B) can easily be positive in one cell but negative in its
    //   coupled cell
    // - so cyclic value might evaluate to 0
    // - which upsets the division
    // - note that instead of sign() any other field might have the same
    //   positive and negative coupled values - it is just a lot less likely
    // - problem is that overridden patch types do not propagate - use in-place
    //   operations only

    volScalarField denom
    (
        IOobject
        (
            "denom",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        sign(B)
    );

    // Override interpolated value on interpolated coupled patches
    for (auto& pfld : denom.boundaryFieldRef())
    {
        if
        (
            isA<coupledFvPatchField<scalar>>(pfld)
        && !isA<processorFvPatchField<scalar>>(pfld)
        )
        {
            pfld = pfld.patchInternalField();
        }
    }

    denom *= kappa_*max(mag(B), dimensionedScalar(B.dimensions(), VSMALL));

    // (O:Eq. 26)
    *result1 = // ObukhovLength
        -min
        (
            dimensionedScalar(dimLength, ROOTVGREAT), // neutral stratification
            pow3(*result2)/denom
        );

    return isNew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ObukhovLength::ObukhovLength
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_("U"),
    resultName1_("ObukhovLength"),
    resultName2_("Ustar"),
    rhoRef_(1.0),
    kappa_(0.40),
    beta_
    (
        dimensionedScalar
        (
            dimless/dimTemperature,
            dict.getCheckOrDefault<scalar>
            (
                "beta",
                3e-3,
                [&](const scalar x){ return x > SMALL; }
            )
        )
    ),
    g_
    (
        "g",
        dimLength/sqr(dimTime),
        meshObjects::gravity::New(mesh_.time()).value()
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ObukhovLength::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict),

    UName_ = dict.getOrDefault<word>("U", "U");
    resultName1_ = dict.getOrDefault<word>("ObukhovLength", "ObukhovLength");
    resultName2_ = dict.getOrDefault<word>("Ustar", "Ustar");

    if (UName_ != "U" && resultName1_ == "ObukhovLength")
    {
        resultName1_ += '(' + UName_ + ')';
    }

    if (UName_ != "U" && resultName1_ == "Ustar")
    {
        resultName2_ += '(' + UName_ + ')';
    }

    rhoRef_ = dict.getOrDefault<scalar>("rhoRef", 1.0);
    kappa_ = dict.getOrDefault<scalar>("kappa", 0.4);
    beta_.value() = dict.getOrDefault<scalar>("beta", 3e-3);

    return true;
}


bool Foam::functionObjects::ObukhovLength::execute()
{
    Log << type() << " " << name() << " execute:" << endl;

    bool isNew = false;

    isNew = calcOL();

    if (isNew) Log << " (new)" << nl << endl;

    return true;
}


bool Foam::functionObjects::ObukhovLength::write()
{
    const auto* ioptr1 = mesh_.cfindObject<regIOobject>(resultName1_);
    const auto* ioptr2 = mesh_.cfindObject<regIOobject>(resultName2_);

    if (ioptr1)
    {
        Log << type() << " " << name() << " write:" << nl
            << "    writing field " << ioptr1->name() << nl
            << "    writing field " << ioptr2->name() << endl;

        ioptr1->write();
        ioptr2->write();
    }

    return true;
}


void Foam::functionObjects::ObukhovLength::removeObukhovLength()
{
    mesh_.thisDb().checkOut(resultName1_);
    mesh_.thisDb().checkOut(resultName2_);
}


void Foam::functionObjects::ObukhovLength::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        removeObukhovLength();
    }
}


void Foam::functionObjects::ObukhovLength::movePoints(const polyMesh& m)
{
    if (&m == &mesh_)
    {
        removeObukhovLength();
    }
}


// ************************************************************************* //
