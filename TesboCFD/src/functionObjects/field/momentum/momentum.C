/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "momentum.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSet.H"
#include "cylindricalRotation.H"
#include "emptyPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(momentum, 0);
    addToRunTimeSelectionTable(functionObject, momentum, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::functionObjects::momentum::purgeFields()
{
    obr_.checkOut(scopedName("momentum"));
    obr_.checkOut(scopedName("angularMomentum"));
    obr_.checkOut(scopedName("angularVelocity"));
}


template<class GeoField>
Foam::autoPtr<GeoField>
Foam::functionObjects::momentum::newField
(
    const word& baseName,
    const dimensionSet& dims,
    bool registerObject
) const
{
    return
        autoPtr<GeoField>::New
        (
            IOobject
            (
                scopedName(baseName),
                time_.timeName(),
                mesh_.thisDb(),
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                registerObject
            ),
            mesh_,
            Foam::zero{},  // value
            dims
        );
}


void Foam::functionObjects::momentum::calc()
{
    initialise();

    // Ensure volRegion is properly up-to-date.
    // Purge old fields if anything has changed (eg, mesh size etc)
    if (volRegion::update())
    {
        purgeFields();
    }

    // When field writing is not enabled we need our local storage
    // for the momentum and angular velocity fields
    autoPtr<volVectorField> tmomentum, tAngularMom, tAngularVel;


    // The base fields required
    const auto& U = lookupObject<volVectorField>(UName_);
    const auto* rhoPtr = findObject<volScalarField>(rhoName_);

    const dimensionedScalar rhoRef("rho", dimDensity, rhoRef_);

    // For quantities such as the mass-averaged angular velocity,
    // we would calculate the mass per-cell here.

    // tmp<volScalarField::Internal> tmass =
    // (
    //     rhoPtr
    //   ? (mesh_.V() * (*rhoPtr))
    //   : (mesh_.V() * rhoRef)
    // );


    // Linear momentum
    // ~~~~~~~~~~~~~~~

    auto* pmomentum = getObjectPtr<volVectorField>(scopedName("momentum"));

    if (!pmomentum)
    {
        tmomentum = newField<volVectorField>("momentum", dimVelocity*dimMass);
        pmomentum = tmomentum.get();  // get(), not release()
    }
    auto& momentum = *pmomentum;

    if (rhoPtr)
    {
        momentum.ref() = (U * mesh_.V() * (*rhoPtr));
    }
    else
    {
        momentum.ref() = (U * mesh_.V() * rhoRef);
    }
    momentum.correctBoundaryConditions();


    // Angular momentum
    // ~~~~~~~~~~~~~~~~

    auto* pAngularMom =
        getObjectPtr<volVectorField>(scopedName("angularMomentum"));

    if (hasCsys_ && !pAngularMom)
    {
        tAngularMom =
            newField<volVectorField>("angularMomentum", dimVelocity*dimMass);
        pAngularMom = tAngularMom.get();  // get(), not release()
    }
    else if (!pAngularMom)
    {
        // Angular momentum not requested, but alias to normal momentum
        // to simplify logic when calculating the summations
        pAngularMom = pmomentum;
    }
    auto& angularMom = *pAngularMom;


    // Angular velocity
    // ~~~~~~~~~~~~~~~~

    auto* pAngularVel =
        getObjectPtr<volVectorField>(scopedName("angularVelocity"));

    if (hasCsys_)
    {
        if (!pAngularVel)
        {
            tAngularVel =
                newField<volVectorField>("angularVelocity", dimVelocity);
            pAngularVel = tAngularVel.get();  // get(), not release()
        }
        auto& angularVel = *pAngularVel;


        // Global to local

        angularVel.primitiveFieldRef() =
            csys_.invTransform(mesh_.cellCentres(), U.internalField());

        angularVel.correctBoundaryConditions();

        if (rhoPtr)
        {
            angularMom.ref() = (angularVel * mesh_.V() * (*rhoPtr));
        }
        else
        {
            angularMom.ref() = (angularVel * mesh_.V() * rhoRef);
        }

        angularMom.correctBoundaryConditions();
    }


    // Integrate the selection

    sumMomentum_ = Zero;
    sumAngularMom_ = Zero;

    if (volRegion::useAllCells())
    {
        for (label celli=0; celli < mesh_.nCells(); ++celli)
        {
            sumMomentum_ += momentum[celli];
            sumAngularMom_ += angularMom[celli];
        }
    }
    else
    {
        for (const label celli : cellIDs())
        {
            sumMomentum_ += momentum[celli];
            sumAngularMom_ += angularMom[celli];
        }
    }

    reduce(sumMomentum_, sumOp<vector>());
    reduce(sumAngularMom_, sumOp<vector>());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::momentum::writeFileHeader(Ostream& os)
{
    if (!writeToFile() || writtenHeader_)
    {
        return;
    }

    if (hasCsys_)
    {
        writeHeader(os, "Momentum, Angular Momentum");
        writeHeaderValue(os, "origin", csys_.origin());
        writeHeaderValue(os, "axis", csys_.e3());
    }
    else
    {
        writeHeader(os, "Momentum");
    }

    if (!volRegion::useAllCells())
    {
        writeHeader
        (
            os,
            "Selection " + regionTypeNames_[regionType_]
          + " = " + regionName_
        );
    }

    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(momentum_x momentum_y momentum_z)");

    if (hasCsys_)
    {
        writeTabbed(os, "(momentum_r momentum_rtheta momentum_axis)");
    }

    writeTabbed(os, "volume");
    os  << endl;

    writtenHeader_ = true;
}


void Foam::functionObjects::momentum::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (!foundObject<volVectorField>(UName_))
    {
        FatalErrorInFunction
            << "Could not find U: " << UName_ << " in database"
            << exit(FatalError);
    }


    const auto* pPtr = cfindObject<volScalarField>(pName_);

    if (pPtr && pPtr->dimensions() == dimPressure)
    {
        // Compressible - rho is mandatory

        if (!foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::momentum::writeValues(Ostream& os)
{
    if (log)
    {
        Info<< type() << " " << name() << " write:" << nl;

        Info<< "    Sum of Momentum";

        if (!volRegion::useAllCells())
        {
            Info<< ' ' << regionTypeNames_[regionType_]
                << ' ' << regionName_;
        }

        Info<< " (volume " << volRegion::V() << ')' << nl
            << "        linear  : " << sumMomentum_ << nl;

        if (hasCsys_)
        {
            Info<< "        angular : " << sumAngularMom_ << nl;
        }

        Info<< endl;
    }

    if (writeToFile())
    {
        writeCurrentTime(os);

        os << tab << sumMomentum_;

        if (hasCsys_)
        {
            os << tab << sumAngularMom_;
        }

        os << tab << volRegion::V() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::momentum::momentum
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(mesh_, name, typeName, dict),
    sumMomentum_(Zero),
    sumAngularMom_(Zero),
    UName_(),
    pName_(),
    rhoName_(),
    rhoRef_(1.0),
    csys_(),
    hasCsys_(false),
    writeMomentum_(false),
    writeVelocity_(false),
    writePosition_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::momentum::momentum
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(mesh_, name, typeName, dict),
    sumMomentum_(Zero),
    sumAngularMom_(Zero),
    UName_(),
    pName_(),
    rhoName_(),
    rhoRef_(1.0),
    csys_(),
    hasCsys_(false),
    writeMomentum_(false),
    writeVelocity_(false),
    writePosition_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::momentum::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    volRegion::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    // Optional entries U and p
    UName_ = dict.getOrDefault<word>("U", "U");
    pName_ = dict.getOrDefault<word>("p", "p");
    rhoName_ = dict.getOrDefault<word>("rho", "rho");
    rhoRef_ = dict.getOrDefault<scalar>("rhoRef", 1.0);
    hasCsys_ = dict.getOrDefault("cylindrical", false);

    if (hasCsys_)
    {
        csys_ = coordSystem::cylindrical(dict);
    }

    writeMomentum_ = dict.getOrDefault("writeMomentum", false);
    writeVelocity_ = dict.getOrDefault("writeVelocity", false);
    writePosition_ = dict.getOrDefault("writePosition", false);

    Info<<"Integrating for selection: "
        << regionTypeNames_[regionType_]
        << " (" << regionName_ << ")" << nl;

    if (writeMomentum_)
    {
        Info<< "    Momentum fields will be written" << endl;

        regIOobject::store
        (
            newField<volVectorField>("momentum", dimVelocity*dimMass)
        );

        if (hasCsys_)
        {
            regIOobject::store
            (
                newField<volVectorField>("angularMomentum", dimVelocity*dimMass)
            );
        }
    }

    if (hasCsys_)
    {
        if (writeVelocity_)
        {
            Info<< "    Angular velocity will be written" << endl;

            regIOobject::store
            (
                newField<volVectorField>("angularVelocity", dimVelocity)
            );
        }

        if (writePosition_)
        {
            Info<< "    Angular position will be written" << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::momentum::execute()
{
    calc();

    if (Pstream::master())
    {
        writeFileHeader(file());

        writeValues(file());

        Log << endl;
    }

    // Write state/results information
    setResult("momentum_x", sumMomentum_[0]);
    setResult("momentum_y", sumMomentum_[1]);
    setResult("momentum_z", sumMomentum_[2]);

    setResult("momentum_r", sumAngularMom_[0]);
    setResult("momentum_rtheta", sumAngularMom_[1]);
    setResult("momentum_axis", sumAngularMom_[2]);

    return true;
}


bool Foam::functionObjects::momentum::write()
{
    if (writeMomentum_ || (hasCsys_ && (writeVelocity_ || writePosition_)))
    {
        Log << "Writing fields" << nl;

        const volVectorField* fieldPtr;

        fieldPtr = findObject<volVectorField>(scopedName("momentum"));
        if (fieldPtr) fieldPtr->write();

        fieldPtr = findObject<volVectorField>(scopedName("angularMomentum"));
        if (fieldPtr) fieldPtr->write();

        fieldPtr = findObject<volVectorField>(scopedName("angularVelocity"));
        if (fieldPtr) fieldPtr->write();

        if (hasCsys_ && writePosition_)
        {
            // Clunky, but currently no simple means of handling
            // component-wise conversion and output

            auto cyl_r = newField<volScalarField>("cyl_r", dimLength);
            auto cyl_t = newField<volScalarField>("cyl_theta", dimless);
            auto cyl_z = newField<volScalarField>("cyl_z", dimLength);

            // Internal
            {
                const auto& pts = mesh_.cellCentres();
                const label len = pts.size();

                UList<scalar>& r = cyl_r->primitiveFieldRef(false);
                UList<scalar>& t = cyl_t->primitiveFieldRef(false);
                UList<scalar>& z = cyl_z->primitiveFieldRef(false);

                for (label i=0; i < len; ++i)
                {
                    point p(csys_.localPosition(pts[i]));

                    r[i] = p.x();
                    t[i] = p.y();
                    z[i] = p.z();
                }
            }

            // Boundary
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

            forAll(pbm, patchi)
            {
                if (isA<emptyPolyPatch>(pbm[patchi]))
                {
                    continue;
                }

                const auto& pts = pbm[patchi].faceCentres();
                const label len = pts.size();

                UList<scalar>& r = cyl_r->boundaryFieldRef(false)[patchi];
                UList<scalar>& t = cyl_t->boundaryFieldRef(false)[patchi];
                UList<scalar>& z = cyl_z->boundaryFieldRef(false)[patchi];

                for (label i=0; i < len; ++i)
                {
                    point p(csys_.localPosition(pts[i]));

                    r[i] = p.x();
                    t[i] = p.y();
                    z[i] = p.z();
                }
            }

            cyl_r->write();
            cyl_t->write();
            cyl_z->write();
        }
    }

    return true;
}


void Foam::functionObjects::momentum::updateMesh(const mapPolyMesh& mpm)
{
    volRegion::updateMesh(mpm);
    purgeFields();  // Mesh motion makes calculated fields dubious
}


void Foam::functionObjects::momentum::movePoints(const polyMesh& pm)
{
    volRegion::movePoints(pm);
    purgeFields();  // Mesh motion makes calculated fields dubious
}


// ************************************************************************* //
