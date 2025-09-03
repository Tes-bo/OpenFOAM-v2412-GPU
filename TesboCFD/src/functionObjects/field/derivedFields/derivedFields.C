/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "derivedFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(derivedFields, 0);
    addToRunTimeSelectionTable(functionObject, derivedFields, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::derivedFields::derivedType
>
Foam::functionObjects::derivedFields::knownNames
({
    { derivedType::NONE , "none" },
    { derivedType::MASS_FLUX , "rhoU" },
    { derivedType::TOTAL_PRESSURE , "pTotal" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static bool calc_rhoU
(
    const fvMesh& mesh,
    const word& derivedName,
    const scalar rhoRef
)
{
    // rhoU = rho * U

    const auto* rhoPtr = mesh.findObject<volScalarField>("rho");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    volVectorField* result = mesh.getObjectPtr<volVectorField>(derivedName);

    const bool isNew = !result;

    if (!result)
    {
        result = new volVectorField
        (
            IOobject
            (
                derivedName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            (dimDensity * dimVelocity)
        );

        result->store();
    }

    if (rhoPtr)
    {
        const auto& rho = *rhoPtr;

        *result = (rho * U);
    }
    else
    {
        const dimensionedScalar rho("rho", dimDensity, rhoRef);

        *result = (rho * U);
    }

    return isNew;
}


static bool calc_pTotal
(
    const fvMesh& mesh,
    const word& derivedName,
    const scalar rhoRef
)
{
    // pTotal = p + rho * U^2 / 2

    const auto* rhoPtr = mesh.findObject<volScalarField>("rho");
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    volScalarField* result = mesh.getObjectPtr<volScalarField>(derivedName);

    const bool isNew = !result;

    if (!result)
    {
        result = new volScalarField
        (
            IOobject
            (
                derivedName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            dimPressure
        );

        result->store();
    }

    if (rhoPtr)
    {
        const auto& rho = *rhoPtr;

        *result = (p + 0.5 * rho * magSqr(U));
    }
    else
    {
        const dimensionedScalar rho("rho", dimDensity, rhoRef);

        *result = (rho * (p + 0.5 * magSqr(U)));
    }

    return isNew;
}
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::derivedFields::derivedFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    derivedTypes_(),
    rhoRef_(1.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::derivedFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict),

    rhoRef_ = dict.getOrDefault<scalar>("rhoRef", 1);

    wordList derivedNames(dict.get<wordList>("derived"));

    derivedTypes_.resize(derivedNames.size());

    label nbad = 0, ngood = 0;

    for (const word& key : derivedNames)
    {
        derivedTypes_[ngood] = knownNames.lookup(key, derivedType::UNKNOWN);

        switch (derivedTypes_[ngood])
        {
            case derivedType::NONE:
                break;

            case derivedType::UNKNOWN:
            {
                derivedNames[nbad++] = key;
                break;
            }

            default:
            {
                ++ngood;
                break;
            }
        }
    }

    if (nbad)
    {
        WarningInFunction
            << "Ignoring unknown derived names: "
            << SubList<word>(derivedNames, nbad) << nl;
    }

    derivedTypes_.resize(ngood);

    // Output the good names
    forAll(derivedTypes_, i)
    {
        derivedNames[i] = knownNames[derivedTypes_[i]];
    }

    Info<< type() << " derived: "
        << flatOutput(SubList<word>(derivedNames, ngood)) << nl;

    return true;
}


bool Foam::functionObjects::derivedFields::execute()
{
    Log << type() << " calculating:";

    for (const derivedType category : derivedTypes_)
    {
        bool isNew = false;

        switch (category)
        {
            case derivedType::MASS_FLUX:
            {
                isNew = calc_rhoU(mesh_, knownNames[category], rhoRef_);

                Log << "  " << knownNames[category];
                if (isNew) Log << " (new)";
                break;
            }

            case derivedType::TOTAL_PRESSURE:
            {
                isNew = calc_pTotal(mesh_, knownNames[category], rhoRef_);

                Log << "  " << knownNames[category];
                if (isNew) Log << " (new)";
                break;
            }

            default:
                break;
        }
    }

    Log << nl << endl;

    return true;
}


bool Foam::functionObjects::derivedFields::write()
{
    for (const derivedType category : derivedTypes_)
    {
        switch (category)
        {
            case derivedType::NONE:
            case derivedType::UNKNOWN:
                break;

            default:
            {
                const auto* ioptr =
                    mesh_.cfindObject<regIOobject>(knownNames[category]);

                if (ioptr)
                {
                    Log << type() << " " << name() << " write:" << nl
                        << "    writing field " << ioptr->name() << endl;

                    ioptr->write();
                }
                break;
            }
        }
    }

    return true;
}


void Foam::functionObjects::derivedFields::removeDerivedFields()
{
    for (const derivedType category : derivedTypes_)
    {
        mesh_.thisDb().checkOut(knownNames[category]);
    }
}


void Foam::functionObjects::derivedFields::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        removeDerivedFields();
    }
}


void Foam::functionObjects::derivedFields::movePoints(const polyMesh& m)
{
    if (&m == &mesh_)
    {
        removeDerivedFields();
    }
}


// ************************************************************************* //
