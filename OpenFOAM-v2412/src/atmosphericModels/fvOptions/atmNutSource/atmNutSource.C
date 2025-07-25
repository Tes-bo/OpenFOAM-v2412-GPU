/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "atmNutSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmNutSource, 0);
    addToRunTimeSelectionTable(option, atmNutSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmNutSource::atmNutSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    artNutName_(dict.getOrDefault<word>("nut", "artNut")),
    artNut_
    (
        IOobject
        (
            artNutName_,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimLength)/dimTime, Zero)
    )
{
    if (!(artNut_.headerOk()))
    {
        FatalErrorInFunction
            << "Unable to find artificial turbulent viscosity field." << nl
            << "atmNutSource requires an artificial nut field."
            << abort(FatalError);
    }

    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.resize(1);

    const tmp<volScalarField>& tnut = turbPtr->nut();

    if (tnut.is_reference())
    {
        fieldNames_[0] = tnut().name();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find nut field." << nl
            << "atmNutSource requires nut field."
            << abort(FatalError);
    }

    fv::option::resetApplied();

    Log << "    Applying atmNutSource to: " << fieldNames_[0] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmNutSource::correct(volScalarField& field)
{
    Log << this->name() << ": correcting " << field.name() << endl;

    field += artNut_;

    field.correctBoundaryConditions();
}


// ************************************************************************* //
