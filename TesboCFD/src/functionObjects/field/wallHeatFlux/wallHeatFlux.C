/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "wallHeatFlux.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "multiphaseInterSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::wallHeatFlux::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& wallHeatFlux
)
{
    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    for (const label patchi : patchIDs_)
    {
        wallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
    }


    const auto* qrPtr = cfindObject<volScalarField>(qrName_);

    if (qrPtr)
    {
        const volScalarField::Boundary& radHeatFluxBf = qrPtr->boundaryField();

        for (const label patchi : patchIDs_)
        {
            wallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFlux::wallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    qrName_("qr")
{
    read(dict);

    writeFileHeader(file());

    volScalarField* wallHeatFluxPtr
    (
        new volScalarField
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
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    regIOobject::store(wallHeatFluxPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFlux::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readIfPresent("qr", qrName_);

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
                << "Requested wall heat-flux on ("
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


bool Foam::functionObjects::wallHeatFlux::execute()
{
    auto& wallHeatFlux = lookupObjectRef<volScalarField>(scopedName(typeName));

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.transport().he(),
            wallHeatFlux
        );
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.he(),
            wallHeatFlux
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.he(), wallHeatFlux);
    }
    else if
    (
        foundObject<multiphaseInterSystem>
            (multiphaseInterSystem::phasePropertiesName)
    )
    {
        const auto& thermo = lookupObject<multiphaseInterSystem>
        (
            multiphaseInterSystem::phasePropertiesName
        );

        calcHeatFlux(thermo.kappaEff()(), thermo.T(), wallHeatFlux);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }


    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();

    for (const label patchi : patchIDs_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFlux.boundaryField()[patchi];

        const MinMax<scalar> limits = gMinMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << limits.min()
                << token::TAB << limits.max()
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << limits.min() << ", " << limits.max()
            << ", " << integralHfp << endl;

        this->setResult("min(" + pp.name() + ")", limits.min());
        this->setResult("max(" + pp.name() + ")", limits.max());
        this->setResult("int(" + pp.name() + ")", integralHfp);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFlux::write()
{
    const auto& wallHeatFlux =
        lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << ' ' << name() << " write:" << nl
        << "    writing field " << wallHeatFlux.name() << endl;

    wallHeatFlux.write();

    return true;
}


// ************************************************************************* //
