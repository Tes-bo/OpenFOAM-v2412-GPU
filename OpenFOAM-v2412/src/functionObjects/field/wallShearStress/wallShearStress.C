/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "wallShearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallShearStress, 0);
    addToRunTimeSelectionTable(functionObject, wallShearStress, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallShearStress::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall shear stress");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os  << endl;
}


void Foam::functionObjects::wallShearStress::calcShearStress
(
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    for (const label patchi : patchIDs_)
    {
        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStress::wallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    writeFields_(true)  // May change in the future
{
    read(dict);

    writeFileHeader(file());

    volVectorField* wallShearStressPtr
    (
        new volVectorField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            mesh_,
            dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    regIOobject::store(wallShearStressPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallShearStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    writeFields_ = true;   // May change in the future
    dict.readIfPresent("writeFields", writeFields_);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    wordRes patchNames;
    labelHashSet patchSet;
    if (dict.readIfPresent("patches", patchNames) && !patchNames.empty())
    {
        patchSet = pbm.patchSet(patchNames);
    }

    labelHashSet allWalls(pbm.findPatchIDs<wallPolyPatch>());

    Info<< type() << ' ' << name() << ':' << nl;

    if (patchSet.empty())
    {
        patchIDs_ = allWalls.sortedToc();

        Info<< "    processing all (" << patchIDs_.size()
            << ") wall patches" << nl << endl;
    }
    else
    {
        allWalls &= patchSet;
        patchSet -= allWalls;
        patchIDs_ = allWalls.sortedToc();

        if (!patchSet.empty())
        {
            WarningInFunction
                << "Requested wall shear stress on ("
                << patchSet.size() << ") non-wall patches:" << nl;

            for (const label patchi : patchSet.sortedToc())
            {
                Info<< "        " << pbm[patchi].name() << nl;
            }
            Info<< nl;
        }

        Info<< "    processing (" << patchIDs_.size()
            << ") wall patches:" << nl;

        for (const label patchi : patchIDs_)
        {
            Info<< "        " << pbm[patchi].name() << nl;
        }
        Info<< endl;
    }

    return true;
}


bool Foam::functionObjects::wallShearStress::execute()
{
    auto& wallShearStress =
        mesh_.lookupObjectRef<volVectorField>(scopedName(typeName));

    // Compressible
    {
        typedef compressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devRhoReff(), wallShearStress);
            return true;
        }
    }

    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devReff(), wallShearStress);
            return true;
        }
    }

    FatalErrorInFunction
        << "Unable to find turbulence model in the "
        << "database" << exit(FatalError);

    return false;
}


bool Foam::functionObjects::wallShearStress::write()
{
    const auto& wallShearStress =
        obr_.lookupObject<volVectorField>(scopedName(typeName));

    Log << type() << ' ' << name() << " write:" << nl;

    if (writeFields_)
    {
        Log << "    writing field " << wallShearStress.name() << endl;
        wallShearStress.write();
    }

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchIDs_)
    {
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = wallShearStress.boundaryField()[patchi];

        const MinMax<vector> limits = gMinMax(ssp);

        if (UPstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << limits.min()
                << token::TAB << limits.max()
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << limits.min() << ", " << limits.max() << endl;
    }

    return true;
}


// ************************************************************************* //
