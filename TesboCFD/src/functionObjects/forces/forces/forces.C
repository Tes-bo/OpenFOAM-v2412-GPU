/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "forces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);
    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forces::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    point origin(Zero);

    // With objectRegistry for access to indirect (global) coordinate systems
    coordSysPtr_ = coordinateSystem::NewIfPresent(obr_, dict);

    if (coordSysPtr_)
    {
        // Report ...
    }
    else if (dict.readIfPresent("CofR", origin))
    {
        const vector e3
        (
            e3Name.empty() ? vector(0, 0, 1) : dict.get<vector>(e3Name)
        );
        const vector e1
        (
            e1Name.empty() ? vector(1, 0, 0) : dict.get<vector>(e1Name)
        );

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        // No 'coordinateSystem' or 'CofR'
        // - enforce a cartesian system

        coordSysPtr_.reset(new coordSystem::cartesian(dict));
    }
}


Foam::volVectorField& Foam::functionObjects::forces::force()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("force"));

    if (!ptr)
    {
        ptr = new volVectorField
        (
            IOobject
            (
                scopedName("force"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimForce, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


Foam::volVectorField& Foam::functionObjects::forces::moment()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("moment"));

    if (!ptr)
    {
        ptr = new volVectorField
        (
            IOobject
            (
                scopedName("moment"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimForce*dimLength, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


void Foam::functionObjects::forces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)
        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_
                << " or p:" << pName_ << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_ << " in database"
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::forces::reset()
{
    sumPatchForcesP_ = Zero;
    sumPatchForcesV_ = Zero;
    sumPatchMomentsP_ = Zero;
    sumPatchMomentsV_ = Zero;

    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    auto& force = this->force();
    auto& moment = this->moment();

    if (porosity_)
    {
        force == dimensionedVector(force.dimensions(), Zero);
        moment == dimensionedVector(moment.dimensions(), Zero);
    }
    else
    {
        constexpr bool updateAccessTime = false;
        for (const label patchi : patchIDs_)
        {
            force.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
            moment.boundaryFieldRef(updateAccessTime)[patchi] = Zero;
        }
    }
}


Foam::tmp<Foam::symmTensorField>
Foam::functionObjects::forces::devRhoReff
(
    const tensorField& gradUp,
    const label patchi
) const
{
    typedef incompressible::turbulenceModel icoTurbModel;
    typedef compressible::turbulenceModel cmpTurbModel;

    if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return -rho(patchi)*turb.nuEff(patchi)*devTwoSymm(gradUp);
    }
    else if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return -turb.muEff(patchi)*devTwoSymm(gradUp);
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(fluidThermo::dictName);

        return -thermo.mu(patchi)*devTwoSymm(gradUp);
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return -rho(patchi)*laminarT.nu(patchi)*devTwoSymm(gradUp);
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return -rho(patchi)*nu.value()*devTwoSymm(gradUp);
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
             lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return volScalarField::New
        (
            "rho",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return (lookupObject<volScalarField>(rhoName_));
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::forces::rho(const label patchi) const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<scalarField>::New
        (
            mesh_.boundary()[patchi].size(),
            rhoRef_
        );
    }

    const auto& rho = lookupObject<volScalarField>(rhoName_);
    return rho.boundaryField()[patchi];
}


Foam::scalar Foam::functionObjects::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::forces::addToPatchFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fP,
    const vectorField& fV
)
{
    constexpr bool updateAccessTime = false;

    sumPatchForcesP_ += sum(fP);
    sumPatchForcesV_ += sum(fV);
    force().boundaryFieldRef(updateAccessTime)[patchi] += fP + fV;

    const vectorField mP(Md^fP);
    const vectorField mV(Md^fV);

    sumPatchMomentsP_ += sum(mP);
    sumPatchMomentsV_ += sum(mV);
    moment().boundaryFieldRef(updateAccessTime)[patchi] += mP + mV;
}


void Foam::functionObjects::forces::addToInternalField
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& f
)
{
    auto& force = this->force();
    auto& moment = this->moment();

    forAll(cellIDs, i)
    {
        const label celli = cellIDs[i];

        sumInternalForces_ += f[i];
        force[celli] += f[i];

        const vector m(Md[i]^f[i]);
        sumInternalMoments_ += m;
        moment[celli] = m;
    }
}


void Foam::functionObjects::forces::createIntegratedDataFiles()
{
    if (!forceFilePtr_)
    {
        forceFilePtr_ = newFileAtStartTime("force");
        writeIntegratedDataFileHeader("Force", forceFilePtr_());
    }

    if (!momentFilePtr_)
    {
        momentFilePtr_ = newFileAtStartTime("moment");
        writeIntegratedDataFileHeader("Moment", momentFilePtr_());
    }
}


void Foam::functionObjects::forces::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    const auto vecDesc = [](const word& root)->string
    {
        return root + "_x " + root + "_y " + root + "_z";
    };
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, vecDesc("total"));
    writeTabbed(os, vecDesc("pressure"));
    writeTabbed(os, vecDesc("viscous"));

    if (porosity_)
    {
        writeTabbed(os, vecDesc("porous"));
    }

    os  << endl;
}


void Foam::functionObjects::forces::writeIntegratedDataFiles()
{
    const auto& coordSys = coordSysPtr_();

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchForcesP_),
        coordSys.localVector(sumPatchForcesV_),
        coordSys.localVector(sumInternalForces_),
        forceFilePtr_()
    );

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchMomentsP_),
        coordSys.localVector(sumPatchMomentsV_),
        coordSys.localVector(sumInternalMoments_),
        momentFilePtr_()
    );
}


void Foam::functionObjects::forces::writeIntegratedDataFile
(
    const vector& pres,
    const vector& vis,
    const vector& internal,
    OFstream& os
) const
{
    writeCurrentTime(os);

    writeValue(os, pres + vis + internal);
    writeValue(os, pres);
    writeValue(os, vis);

    if (porosity_)
    {
        writeValue(os, internal);
    }

    os  << endl;
}


void Foam::functionObjects::forces::logIntegratedData
(
    const string& descriptor,
    const vector& pres,
    const vector& vis,
    const vector& internal
) const
{
    if (!log)
    {
        return;
    }

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << (pres + vis + internal) << nl
        << "        Pressure : " << pres << nl
        << "        Viscous  : " << vis << nl;

    if (porosity_)
    {
        Log << "        Porous   : " << internal << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << ' ' << name() << ':' << endl;

    // Can also use pbm.indices(), but no warnings...
    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();

    dict.readIfPresent("directForceDensity", directForceDensity_);
    if (directForceDensity_)
    {
        // Optional name entry for fD
        if (dict.readIfPresent<word>("fD", fDName_))
        {
            Info<< "    fD: " << fDName_ << endl;
        }
    }
    else
    {
        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }
        if (dict.readIfPresent<word>("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
            Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent<scalar>("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
        }
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }


    return true;
}


void Foam::functionObjects::forces::calcForcesMoments()
{
    initialise();

    reset();

    const point& origin = coordSysPtr_->origin();

    if (directForceDensity_)
    {
        const auto& fD = lookupObject<volVectorField>(fDName_);
        const auto& fDb = fD.boundaryField();

        const auto& Sfb = mesh_.Sf().boundaryField();
        const auto& magSfb = mesh_.magSf().boundaryField();
        const auto& Cb = mesh_.C().boundaryField();

        for (const label patchi : patchIDs_)
        {
            const vectorField Md(Cb[patchi] - origin);

            // Pressure force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fP
            (
                Sfb[patchi]/magSfb[patchi]
               *(
                    Sfb[patchi] & fDb[patchi]
                )
            );

            // Viscous force (total force minus pressure fP)
            const vectorField fV(magSfb[patchi]*fDb[patchi] - fP);

            addToPatchFields(patchi, Md, fP, fV);
        }
    }
    else
    {
        const auto& p = lookupObject<volScalarField>(pName_);
        const auto& pb = p.boundaryField();

        const auto& Sfb = mesh_.Sf().boundaryField();
        const auto& Cb = mesh_.C().boundaryField();

        const auto& U = lookupObject<volVectorField>(UName_);
        tmp<volTensorField> tgradU = fvc::grad(U);
        const volTensorField& gradU = tgradU();
        const auto& gradUb = gradU.boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar rhoRef = rho(p);
        const scalar pRef = pRef_/rhoRef;

        for (const label patchi : patchIDs_)
        {
            const vectorField Md(Cb[patchi] - origin);

            const vectorField fP(rhoRef*Sfb[patchi]*(pb[patchi] - pRef));

            const vectorField fV
            (
                Sfb[patchi] & devRhoReff(gradUb[patchi], patchi)
            );

            addToPatchFields(patchi, Md, fP, fV);
        }
    }

    if (porosity_)
    {
        const auto& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const UPtrList<const porosityModel> models
        (
            obr_.csorted<porosityModel>()
        );

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, "
                << "but no porosity models found in the database"
                << endl;
        }

        for (const porosityModel& mdl : models)
        {
            // Non-const access required if mesh is changing
            auto& pm = const_cast<porosityModel&>(mdl);

            const vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - origin);

                addToInternalField(cZone, Md, fP);
            }
        }
    }

    reduce(sumPatchForcesP_, sumOp<vector>());
    reduce(sumPatchForcesV_, sumOp<vector>());
    reduce(sumPatchMomentsP_, sumOp<vector>());
    reduce(sumPatchMomentsV_, sumOp<vector>());
    reduce(sumInternalForces_, sumOp<vector>());
    reduce(sumInternalMoments_, sumOp<vector>());
}


Foam::vector Foam::functionObjects::forces::forceEff() const
{
    return sumPatchForcesP_ + sumPatchForcesV_ + sumInternalForces_;
}


Foam::vector Foam::functionObjects::forces::momentEff() const
{
    return sumPatchMomentsP_ + sumPatchMomentsV_ + sumInternalMoments_;
}


bool Foam::functionObjects::forces::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    const auto localFp(coordSys.localVector(sumPatchForcesP_));
    const auto localFv(coordSys.localVector(sumPatchForcesV_));
    const auto localFi(coordSys.localVector(sumInternalForces_));

    logIntegratedData("forces", localFp, localFv, localFi);

    const auto localMp(coordSys.localVector(sumPatchMomentsP_));
    const auto localMv(coordSys.localVector(sumPatchMomentsV_));
    const auto localMi(coordSys.localVector(sumInternalMoments_));

    logIntegratedData("moments", localMp, localMv, localMi);

    setResult("pressureForce", localFp);
    setResult("viscousForce", localFv);
    setResult("internalForce", localFi);
    setResult("pressureMoment", localMp);
    setResult("viscousMoment", localMv);
    setResult("internalMoment", localMi);

    return true;
}


bool Foam::functionObjects::forces::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment files." << endl;

        createIntegratedDataFiles();
        writeIntegratedDataFiles();
    }

    if (writeFields_)
    {
        Log << "    writing force and moment fields." << endl;

        force().write();
        moment().write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
