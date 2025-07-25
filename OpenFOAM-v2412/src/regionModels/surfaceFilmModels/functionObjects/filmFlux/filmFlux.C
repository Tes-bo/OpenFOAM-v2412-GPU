/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "filmFlux.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(filmFlux, 0);
    addToRunTimeSelectionTable(functionObject, filmFlux, dictionary);
}
}


const Foam::functionObjects::filmFlux::filmType&
Foam::functionObjects::filmFlux::filmModel()
{
    return time_.objectRegistry::lookupObject<filmType>(filmName_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::filmFlux::filmFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    filmName_("surfaceFilmProperties"),
    resultName_(scopedName("filmFlux"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::filmFlux::read(const dictionary& dict)
{
    if (stateFunctionObject::read(dict))
    {
        dict.readIfPresent<word>("film", filmName_);
        dict.readIfPresent<word>("result", resultName_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::filmFlux::execute()
{
    const auto& model = filmModel();

    const fvMesh& filmMesh = model.regionMesh();

    auto* resultPtr = filmMesh.getObjectPtr<surfaceScalarField>(resultName_);

    if (!resultPtr)
    {
        resultPtr = new surfaceScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                filmMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            filmMesh,
            dimensionedScalar(dimMass/dimTime, Zero)
        );

        regIOobject::store(resultPtr);
    }

    auto& result = *resultPtr;

    const surfaceScalarField& phi = model.phi();
    const volScalarField& magSf = model.magSf();
    const volScalarField::Internal& vol = filmMesh.V();

    volScalarField height
    (
        IOobject
        (
            "height",
            time_.timeName(),
            filmMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        filmMesh,
        dimensionedScalar(dimLength, Zero),
        fvPatchFieldBase::zeroGradientType()
    );

    auto& heightc = height.internalFieldRef();

    heightc = max(dimensionedScalar("eps", dimLength, ROOTVSMALL), vol/magSf());
    height.correctBoundaryConditions();

    result = phi/fvc::interpolate(height);

    return true;
}


bool Foam::functionObjects::filmFlux::write()
{
    const auto& filmMesh = filmModel().regionMesh();

    filmMesh.lookupObject<surfaceScalarField>(resultName_).write();

    return true;
}


// ************************************************************************* //
