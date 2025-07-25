/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "faOption.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fa
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fa::option::resetApplied()
{
    applied_.resize_nocopy(fieldNames_.size());
    applied_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    modelType_(modelType),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    fieldNames_(),
    applied_(),
    regionName_(dict.get<word>("region")),
    regionMeshPtr_(nullptr),
    vsmPtr_(nullptr),
    active_(dict.getOrDefault("active", true)),
    log(true)
{
    Log << incrIndent << indent << "Source: " << name_ << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fa::option> Foam::fa::option::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    const word modelType(coeffs.get<word>("type"));

    Info<< indent
        << "Selecting finite area options type " << modelType << endl;

    mesh.time().libs().open
    (
        coeffs,
        "libs",
        dictionaryConstructorTablePtr_
    );

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            coeffs,
            "faOption",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<fa::option>(ctorPtr(name, modelType, coeffs, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fa::option::isActive()
{
    return active_;
}


Foam::label Foam::fa::option::applyToField(const word& fieldName) const
{
    return fieldNames_.find(fieldName);
}


void Foam::fa::option::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningInFunction
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    faMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    faMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    faMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    faMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::constrain(faMatrix<scalar>& eqn, const label fieldi)
{}


void Foam::fa::option::constrain(faMatrix<vector>& eqn, const label fieldi)
{}


void Foam::fa::option::constrain
(
    faMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::constrain
(
    faMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fa::option::constrain(faMatrix<tensor>& eqn, const label fieldi)
{}


void Foam::fa::option::correct(areaScalarField& field)
{}


void Foam::fa::option::correct(areaVectorField& field)
{}


void Foam::fa::option::correct(areaSphericalTensorField& field)
{}


void Foam::fa::option::correct(areaSymmTensorField& field)
{}


void Foam::fa::option::correct(areaTensorField& field)
{}


// ************************************************************************* //
