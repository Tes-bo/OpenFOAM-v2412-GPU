/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "momentumError.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(momentumError, 0);
    addToRunTimeSelectionTable(functionObject, momentumError, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::functionObjects::momentumError::divDevRhoReff()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const auto& U = lookupObject<volVectorField>(UName_);
    tmp<volVectorField> tU(U);

    {
        auto* turb = findObject<cmpTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            tmp<volScalarField> trho = turb->rho();
            tmp<volScalarField> tnuEff = turb->nuEff();

            if (zoneSubSetPtr_)
            {
                const fvMeshSubset& subsetter = zoneSubSetPtr_->subsetter();

                tU = subsetter.interpolate(U, false);
                trho = subsetter.interpolate(turb->rho(), false);
                tnuEff = subsetter.interpolate(turb->nuEff()(), false);
            }

            return tmp<volVectorField>::New
            (
                "divDevRhoReff",
                - fvc::div
                (
                    (trho()*tnuEff())
                    *dev2(T(fvc::grad(tU()))),
                    "div(((rho*nuEff)*dev2(T(grad(U)))))"
                )
                - fvc::laplacian
                (
                    trho()*tnuEff(),
                    tU(),
                    "laplacian(nuEff,U)"
                )
            );
        }
    }

    {
        const auto* turb = findObject<icoTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            tmp<volScalarField> tnuEff = turb->nuEff();

            if (zoneSubSetPtr_)
            {
                const fvMeshSubset& subsetter = zoneSubSetPtr_->subsetter();

                tU = subsetter.interpolate(U, false);
                tnuEff = subsetter.interpolate(turb->nuEff()(), false);
            }

            return tmp<volVectorField>::New
            (
                "divDevRhoReff",
                - fvc::div
                (
                    tnuEff()*dev2(T(fvc::grad(tU()))),
                    "div((nuEff*dev2(T(grad(U)))))"
                )
                - fvc::laplacian(tnuEff(), tU(), "laplacian(nuEff,U)")
            );
        }
     }

     return volVectorField::null();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::momentumError::momentumError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_("p"),
    UName_("U"),
    phiName_("phi")
{
    read(dict);

    const auto& phi = lookupObject<surfaceScalarField>(phiName_);

    const dimensionSet momDims
    (
        phi.dimensions()*dimVelocity/dimVolume
    );

    if (zoneSubSetPtr_)
    {
        const fvMesh& subMesh = zoneSubSetPtr_->subsetter().subMesh();

        // momentErrorMap

        volVectorField* fldPtr = new volVectorField
        (
            IOobject
            (
                scopedName("momentErrorMap"),
                subMesh.time().timeName(),
                subMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            subMesh,
            dimensionedVector(momDims, Zero)
        );

        regIOobject::store(fldPtr);

    }

    const word momName =
    (
        zoneSubSetPtr_
      ? scopedName("momentErrorZone")
      : scopedName("momentError")
    );

    volVectorField* fldPtr = new volVectorField
    (
        IOobject
        (
            momName,
            time_.timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedVector(momDims, Zero)
    );

    regIOobject::store(fldPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::momentumError::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << ' ' << name() << ':' << nl;

        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }

        if (dict.readIfPresent<word>("phi", phiName_))
        {
            Info<< "    phi: " << phiName_ << endl;
        }

        if (dict.found("cellZones"))
        {
            zoneSubSetPtr_.reset(new Detail::zoneSubSet(mesh_, dict));
        }
        else
        {
            zoneSubSetPtr_.reset(nullptr);
        }

        return true;
    }

    return false;
}


void Foam::functionObjects::momentumError::calcMomentError()
{
    const auto& p = lookupObject<volScalarField>(pName_);
    const auto& U = lookupObject<volVectorField>(UName_);
    const auto& phi = lookupObject<surfaceScalarField>(phiName_);

    if (zoneSubSetPtr_)
    {
        const fvMeshSubset& subsetter = zoneSubSetPtr_->subsetter();

        fvMesh& subMesh = zoneSubSetPtr_->subsetter().subMesh();

        subMesh.fvSchemes::readOpt(mesh_.fvSchemes::readOpt());
        subMesh.fvSchemes::read();

        auto& momentErrMap =
            subMesh.lookupObjectRef<volVectorField>
            (
                scopedName("momentErrorMap")
            );

        tmp<volScalarField> tp = subsetter.interpolate(p, false);
        tmp<volVectorField> tU = subsetter.interpolate(U, false);
        tmp<surfaceScalarField> tphi = subsetter.interpolate(phi, false);

         momentErrMap =
         (
            divDevRhoReff()
          + fvc::div(tphi, tU, "div(phi,U)")
          + fvc::grad(tp, "grad(p)")
         );
    }
    else
    {
        auto& momentErr =
            lookupObjectRef<volVectorField>(scopedName("momentError"));

        momentErr =  fvc::div(phi, U) + fvc::grad(p) + divDevRhoReff();
    }
}


bool Foam::functionObjects::momentumError::execute()
{
    calcMomentError();

    return true;
}


bool Foam::functionObjects::momentumError::write()
{
    Log << "    functionObjects::" << type() << " " << name();

    if (!zoneSubSetPtr_)
    {
        Log << " writing field: " << scopedName("momentError") << endl;

        const auto& momentErr =
            lookupObjectRef<volVectorField>(scopedName("momentError"));

        momentErr.write();
    }
    else
    {
        Log << " writing field: " << scopedName("momentErrorMap") << endl;

        const fvMeshSubset& subsetter = zoneSubSetPtr_->subsetter();
        const fvMesh& subMesh = subsetter.subMesh();

        const auto& momentErrMap =
            subMesh.lookupObject<volVectorField>
            (
                scopedName("momentErrorMap")
            );

        tmp<volVectorField> mapMomErr =
            zoneSubSetPtr_->mapToZone<vector>(momentErrMap);

        mapMomErr().write();
    }

    return true;
}


// ************************************************************************* //
