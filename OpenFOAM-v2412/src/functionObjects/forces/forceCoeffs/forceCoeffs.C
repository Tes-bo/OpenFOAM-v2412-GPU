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

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::initialise()
{
    if (initialised_)
    {
        return;
    }

    initialised_ = true;
}


Foam::volVectorField& Foam::functionObjects::forceCoeffs::forceCoeff()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("forceCoeff"));

    if (!ptr)
    {
        ptr = new volVectorField
        (
            IOobject
            (
                scopedName("forceCoeff"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimless, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


Foam::volVectorField& Foam::functionObjects::forceCoeffs::momentCoeff()
{
    auto* ptr = mesh_.getObjectPtr<volVectorField>(scopedName("momentCoeff"));

    if (!ptr)
    {
        ptr = new volVectorField
        (
            IOobject
            (
                scopedName("momentCoeff"),
                time_.timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedVector(dimless, Zero)
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


void Foam::functionObjects::forceCoeffs::reset()
{
    Cf_.reset();
    Cm_.reset();

    forceCoeff() == dimensionedVector(dimless, Zero);
    momentCoeff() == dimensionedVector(dimless, Zero);
}


Foam::HashTable<Foam::functionObjects::forceCoeffs::coeffDesc>
Foam::functionObjects::forceCoeffs::selectCoeffs() const
{
    HashTable<coeffDesc> coeffs(16);

    coeffs.insert("Cd", coeffDesc("Drag force", "Cd", 0, 0));
    coeffs.insert("Cs", coeffDesc("Side force", "Cs", 1, 2));
    coeffs.insert("Cl", coeffDesc("Lift force", "Cl", 2, 1));

    // Add front/rear options
    const auto frontRearCoeffs(coeffs);
    forAllConstIters(frontRearCoeffs, iter)
    {
        const auto& fr = iter.val();
        coeffs.insert(fr.frontName(), fr.front());
        coeffs.insert(fr.rearName(), fr.rear());
    }

    // Add moments
    coeffs.insert("CmRoll", coeffDesc("Roll moment", "CmRoll", 0, -1));
    coeffs.insert("CmPitch", coeffDesc("Pitch moment", "CmPitch", 1, -1));
    coeffs.insert("CmYaw", coeffDesc("Yaw moment", "CmYaw", 2, -1));

    return coeffs;
}


void Foam::functionObjects::forceCoeffs::calcForceCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar forceScaling
    (
        dimless/dimForce,
        scalar(1)/(Aref_*pDyn + SMALL)
    );

    const auto& coordSys = coordSysPtr_();

    // Calculate force coefficients
    forceCoeff() = forceScaling*force();

    Cf_.reset
    (
        forceScaling.value()*coordSys.localVector(sumPatchForcesP_),
        forceScaling.value()*coordSys.localVector(sumPatchForcesV_),
        forceScaling.value()*coordSys.localVector(sumInternalForces_)
    );
}


void Foam::functionObjects::forceCoeffs::calcMomentCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar momentScaling
    (
        dimless/(dimForce*dimLength),
        scalar(1)/(Aref_*pDyn*lRef_ + SMALL)
    );

    const auto& coordSys = coordSysPtr_();

    // Calculate moment coefficients
    momentCoeff() = momentScaling*moment();

    Cm_.reset
    (
        momentScaling.value()*coordSys.localVector(sumPatchMomentsP_),
        momentScaling.value()*coordSys.localVector(sumPatchMomentsV_),
        momentScaling.value()*coordSys.localVector(sumInternalMoments_)
    );
}


void Foam::functionObjects::forceCoeffs::createIntegratedDataFile()
{
    if (!coeffFilePtr_)
    {
        coeffFilePtr_ = newFileAtStartTime("coefficient");
        writeIntegratedDataFileHeader("Coefficients", coeffFilePtr_());
    }
}


void Foam::functionObjects::forceCoeffs::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();

    writeHeader(os, "Force and moment coefficients");
    writeHeaderValue(os, "dragDir", coordSys.e1());
    writeHeaderValue(os, "sideDir", coordSys.e2());
    writeHeaderValue(os, "liftDir", coordSys.e3());
    writeHeaderValue(os, "rollAxis", coordSys.e1());
    writeHeaderValue(os, "pitchAxis", coordSys.e2());
    writeHeaderValue(os, "yawAxis", coordSys.e3());
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");

    for (const auto& iter : coeffs_.csorted())
    {
        const auto& coeff = iter.val();

        if (!coeff.active_) continue;

        writeTabbed(os, coeff.name_);
    }

    os  << endl;
}


void Foam::functionObjects::forceCoeffs::writeIntegratedDataFile()
{
    OFstream& os = coeffFilePtr_();

    writeCurrentTime(os);

    for (const auto& iter : coeffs_.csorted())
    {
        const auto& coeff = iter.val();

        if (!coeff.active_) continue;

        const vector coeffValue = coeff.value(Cf_, Cm_);

        os  << tab << (coeffValue.x() + coeffValue.y() + coeffValue.z());
    }

    coeffFilePtr_() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forces(name, runTime, dict, false),
    Cf_(),
    Cm_(),
    coeffFilePtr_(),
    magUInf_(Zero),
    lRef_(Zero),
    Aref_(Zero),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }

    initialised_ = false;

    // Reference scales
    dict.readEntry("magUInf", magUInf_);
    dict.readEntry("lRef", lRef_);
    dict.readEntry("Aref", Aref_);

    // If case is compressible we must read rhoInf
    // (stored in rhoRef_) to calculate the reference dynamic pressure
    // Note: for incompressible, rhoRef_ is already initialised to 1
    if (rhoName_ != "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }

    Info<< "    magUInf: " << magUInf_ << nl
        << "    lRef   : " << lRef_ << nl
        << "    Aref   : " << Aref_ << nl
        << "    rhoInf : " << rhoRef_ << endl;

    if (min(magUInf_, rhoRef_) < SMALL || min(lRef_, Aref_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Non-zero values are required for reference scales"
            << exit(FatalIOError);

        return false;
    }

    if (!dict.found("coefficients"))
    {
        Info<< "    Selecting all coefficients" << nl;

        coeffs_ = selectCoeffs();

        for (auto& iter : coeffs_.sorted())
        {
            auto& coeff = iter.val();
            coeff.active_ = true;
            Info<< "    - " << coeff << nl;
        }
    }
    else
    {
        const wordHashSet coeffs(dict.get<wordHashSet>("coefficients"));

        coeffs_ = selectCoeffs();

        Info<< "    Selecting coefficients:" << nl;

        for (const word& key : coeffs)
        {
            auto iter = coeffs_.find(key);

            if (!iter.good())
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown coefficient type " << key
                    << " . Valid entries are : " << coeffs_.sortedToc()
                    << exit(FatalIOError);
            }
            auto& coeff = iter.val();
            coeff.active_ = true;
            Info<< "    - " << coeff << nl;
        }
    }

    Info<< endl;

    return true;
}



bool Foam::functionObjects::forceCoeffs::execute()
{
    forces::calcForcesMoments();

    initialise();

    reset();

    Log << type() << " " << name() << " write:" << nl
        << "    " << "Coefficient"
        << tab << "Total"
        << tab << "Pressure"
        << tab << "Viscous"
        << tab << "Internal"
        << nl;

    calcForceCoeffs();

    calcMomentCoeffs();

    auto logValues = [](const word& name, const vector& coeff, Ostream& os)
    {
        os  << "    " << name << ":"
            << tab << coeff.x() + coeff.y() + coeff.z()
            << tab << coeff.x()
            << tab << coeff.y()
            << tab << coeff.z()
            << nl;
    };

    // Always setting all results
    for (const auto& iter : coeffs_.csorted())
    {
        const word& coeffName = iter.key();
        const auto& coeff = iter.val();

        // Vectors for x:pressure, y:viscous, z:internal
        const vector coeffValue = coeff.value(Cf_, Cm_);

        if (log && coeff.active_)
        {
            logValues(coeffName, coeffValue, Info);
        }

        setResult(coeffName, coeffValue.x() + coeffValue.y() + coeffValue.z());
        setResult(coeffName & "pressure", coeffValue.x());
        setResult(coeffName & "viscous", coeffValue.y());
        setResult(coeffName & "internal", coeffValue.z());
    }

    Log  << endl;

    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment coefficient files." << endl;

        createIntegratedDataFile();

        writeIntegratedDataFile();
    }

    if (writeFields_)
    {
        forceCoeff().write();
        momentCoeff().write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
