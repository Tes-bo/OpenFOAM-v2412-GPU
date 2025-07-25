/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
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

#include "atmPlantCanopyUSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmPlantCanopyUSource, 0);
    addToRunTimeSelectionTable(option, atmPlantCanopyUSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::volScalarField& Foam::fv::atmPlantCanopyUSource::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = mesh_.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmPlantCanopyUSource::atmPlantCanopyUSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    CdName_(),
    LADname_()
{
    read(dict);

    fieldNames_.resize(1, "U");

    fv::option::resetApplied();

    Log << "    Applying atmPlantCanopyUSource to: " << fieldNames_[0] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmPlantCanopyUSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V_ < VSMALL)
    {
        return;
    }

    const volVectorField& U = eqn.psi();
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    // (SP:Eq. 42), (BSG:Eq. 7)
    eqn -= fvm::Sp(Cd*LAD*mag(U), U);
}


void Foam::fv::atmPlantCanopyUSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V_ < VSMALL)
    {
        return;
    }

    const volVectorField& U = eqn.psi();
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    eqn -= fvm::Sp(rho*Cd*LAD*mag(U), U);
}


void Foam::fv::atmPlantCanopyUSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V_ < VSMALL)
    {
        return;
    }

    const volVectorField& U = eqn.psi();
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    eqn -= fvm::Sp(alpha*rho*Cd*LAD*mag(U), U);
}


bool Foam::fv::atmPlantCanopyUSource::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    CdName_ = dict.getOrDefault<word>("Cd", "Cd");
    LADname_ = dict.getOrDefault<word>("LAD", "LAD");

    (void) getOrReadField(CdName_);
    (void) getOrReadField(LADname_);

    return true;
}


// ************************************************************************* //
