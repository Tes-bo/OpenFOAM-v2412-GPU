/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "multiphaseMangrovesTurbulenceModel.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseMangrovesTurbulenceModel, 0);
    addToRunTimeSelectionTable
    (
        option,
        multiphaseMangrovesTurbulenceModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::multiphaseMangrovesTurbulenceModel::kCoeff
(
    const volVectorField& U
) const
{
    auto tkCoeff = volScalarField::New
    (
        IOobject::scopedName(typeName, "kCoeff"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& kCoeff = tkCoeff.ref();

    forAll(zoneIDs_, i)
    {
        const scalar a = aZone_[i];
        const scalar N = NZone_[i];
        const scalar Ckp = CkpZone_[i];
        const scalar Cd = CdZone_[i];

        for (label zonei : zoneIDs_[i])
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (label celli : cz)
            {
                kCoeff[celli] = Ckp*Cd*a*N*mag(U[celli]);
            }
        }
    }

    kCoeff.correctBoundaryConditions();

    return tkCoeff;
}


Foam::tmp<Foam::volScalarField>
Foam::fv::multiphaseMangrovesTurbulenceModel::epsilonCoeff
(
    const volVectorField& U
) const
{
    auto tepsilonCoeff = volScalarField::New
    (
        IOobject::scopedName(typeName, "epsilonCoeff"),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& epsilonCoeff = tepsilonCoeff.ref();

    forAll(zoneIDs_, i)
    {
        const scalar a = aZone_[i];
        const scalar N = NZone_[i];
        const scalar Cep = CepZone_[i];
        const scalar Cd = CdZone_[i];

        for (label zonei : zoneIDs_[i])
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (label celli : cz)
            {
                epsilonCoeff[celli] = Cep*Cd*a*N*mag(U[celli]);
            }
        }
    }

    epsilonCoeff.correctBoundaryConditions();

    return tepsilonCoeff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseMangrovesTurbulenceModel::multiphaseMangrovesTurbulenceModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    aZone_(),
    NZone_(),
    CkpZone_(),
    CepZone_(),
    CdZone_(),
    UName_("U"),
    kName_("k"),
    epsilonName_("epsilon")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseMangrovesTurbulenceModel::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    if (eqn.psi().name() == epsilonName_)
    {
        fvMatrix<scalar> epsilonEqn
        (
          - fvm::Sp(epsilonCoeff(U), eqn.psi())
        );
        eqn += epsilonEqn;
    }
    else if (eqn.psi().name() == kName_)
    {
        fvMatrix<scalar> kEqn
        (
          - fvm::Sp(kCoeff(U), eqn.psi())
        );
        eqn += kEqn;
    }
}


void Foam::fv::multiphaseMangrovesTurbulenceModel::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    if (eqn.psi().name() == epsilonName_)
    {
        fvMatrix<scalar> epsilonEqn
        (
          - fvm::Sp(rho*epsilonCoeff(U), eqn.psi())
        );
        eqn += epsilonEqn;
    }
    else if (eqn.psi().name() == kName_)
    {
        fvMatrix<scalar> kEqn
        (
          - fvm::Sp(rho*kCoeff(U), eqn.psi())
        );
        eqn += kEqn;
    }
}


bool Foam::fv::multiphaseMangrovesTurbulenceModel::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        if (!coeffs_.readIfPresent("epsilonNames", fieldNames_))
        {
            if (coeffs_.found("epsilon"))
            {
                fieldNames_.resize(1);
                coeffs_.readEntry("epsilon", fieldNames_.first());
            }
            else
            {
                fieldNames_.resize(2);
                fieldNames_[0] = "epsilon";
                fieldNames_[1] = "k";
            }
        }
        fv::option::resetApplied();

        // Create the Mangroves models - 1 per region
        const dictionary& regionsDict(coeffs_.subDict("regions"));
        const wordList regionNames(regionsDict.toc());
        aZone_.setSize(regionNames.size(), 1);
        NZone_.setSize(regionNames.size(), 1);
        CkpZone_.setSize(regionNames.size(), 1);
        CepZone_.setSize(regionNames.size(), 1);
        CdZone_.setSize(regionNames.size(), 1);
        zoneIDs_.setSize(regionNames.size());

        forAll(zoneIDs_, i)
        {
            const word& regionName = regionNames[i];
            const dictionary& modelDict = regionsDict.subDict(regionName);

            const word zoneName(modelDict.get<word>("cellZone"));

            zoneIDs_[i] = mesh_.cellZones().indices(zoneName);
            if (zoneIDs_[i].empty())
            {
                FatalErrorInFunction
                    << "Unable to find cellZone " << zoneName << nl
                    << "Valid cellZones are:" << mesh_.cellZones().names()
                    << exit(FatalError);
            }

            modelDict.readEntry("a", aZone_[i]);
            modelDict.readEntry("N", NZone_[i]);
            modelDict.readEntry("Ckp", CkpZone_[i]);
            modelDict.readEntry("Cep", CepZone_[i]);
            modelDict.readEntry("Cd", CdZone_[i]);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
