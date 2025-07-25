/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "surfaceCourantNumber.H"
#include "fvMesh.H"
#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "facEdgeIntegrate.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceCourantNumber, 0);
    addToRunTimeSelectionTable(functionObject, surfaceCourantNumber, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::surfaceCourantNumber::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Surface Courant Number");

    writeCommented(os, "Time");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "mean");
    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceCourantNumber::surfaceCourantNumber
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    resultName_("surfaceCo"),
    phisName_("phis"),
    rhoName_("rho")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfaceCourantNumber::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    dict.readIfPresent("result", resultName_);
    dict.readIfPresent("phis", phisName_);
    dict.readIfPresent("rho", rhoName_);

    // Registry containing all finite-area meshes on the polyMesh
    const auto* faRegistryPtr = faMesh::registry(mesh_);

    if (!faRegistryPtr)
    {
        FatalIOErrorInFunction(dict)
            << "No finite-area object registry is available."
            << abort(FatalIOError);
    }

    word areaName;

    if (!dict.readIfPresent("area", areaName))
    {
        wordList available = faRegistryPtr->sortedNames<faMesh>();
        if (!available.empty())
        {
            areaName = available.front();
        }
    }

    if (areaName.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No name for finite-area mesh is available."
            << abort(FatalIOError);
    }

    faMeshPtr_ = std::shared_ptr<const faMesh>
    (
        faRegistryPtr->cfindObject<faMesh>(areaName),
        [](const faMesh*) { /* no-op deleter to avoid double deletion */ }
    );

    return true;
}


bool Foam::functionObjects::surfaceCourantNumber::execute()
{
    if (!faMeshPtr_->foundObject<edgeScalarField>(phisName_))
    {
        WarningInFunction
            << "No edge flux field is available. "
            << "Name of provided edge flux field (phi): " << phisName_
            << endl;

        return false;
    }

    const auto& phis = faMeshPtr_->lookupObject<edgeScalarField>(phisName_);

    tmp<areaScalarField::Internal> tCo =
        (0.5*faMeshPtr_->time().deltaT())
       *fac::edgeSum(mag(phis))()()
       /faMeshPtr_->S();

    areaScalarField::Internal Co = tCo.ref();

    if (Co.dimensions() == dimDensity)
    {
        Co /= faMeshPtr_->lookupObject<areaScalarField>(rhoName_);
    }

    auto* resultPtr = faMeshPtr_->getObjectPtr<areaScalarField>(resultName_);

    if (!resultPtr)
    {
        resultPtr = new areaScalarField
        (
            IOobject
            (
                resultName_,
                faMeshPtr_->time().timeName(),
                *faMeshPtr_,
                IOobjectOption()
            ),
            *faMeshPtr_,
            dimensionedScalar(dimless, Zero),
            faPatchFieldBase::zeroGradientType()
        );
        regIOobject::store(resultPtr);
    }
    auto& result = *resultPtr;

    result.internalFieldRef() = tCo;
    result.correctBoundaryConditions();


    const scalarMinMax limits(gMinMax(result));
    const scalar mean = gAverage(result);

    Log << "Surface Courant number: "
        << "mean: " << mean
        << " max: " << limits.max()
        << endl;

    if (writeToFile())
    {
        if (!writtenHeader_) writeFileHeader(file());

        writeCurrentTime(file());
        file()
            << token::TAB << limits.min()
            << token::TAB << limits.max()
            << token::TAB << mean
            << endl;
    }

    return true;
}


bool Foam::functionObjects::surfaceCourantNumber::write()
{
    const auto* result = faMeshPtr_->cfindObject<areaScalarField>(resultName_);

    if (!result)
    {
        return false;
    }

    Log << type() << " " << name() << " write: " << result->name() << endl;

    result->write();

    return true;
}


// ************************************************************************* //
