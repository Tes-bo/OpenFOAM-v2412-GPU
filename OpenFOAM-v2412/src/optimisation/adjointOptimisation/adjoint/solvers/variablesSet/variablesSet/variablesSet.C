/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "variablesSet.H"
#include "surfaceFields.H"
#include "fixedValueFvPatchFields.H"
#include "linear.H"

#include "wallFvPatch.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "processorFvPatch.H"
#include "processorFvPatchField.H"
#include "cyclicFvPatch.H"
#include "cyclicFvPatchField.H"
#include "cyclicAMIFvPatch.H"
#include "cyclicAMIFvPatchField.H"
#include "symmetryFvPatch.H"
#include "symmetryFvPatchField.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryPlaneFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(variablesSet, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

variablesSet::variablesSet
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    solverName_(dict.dictName()),
    useSolverNameForFields_
    (
        dict.getOrDefault<bool>("useSolverNameForFields", false)
    )
{}


autoPtr<variablesSet> variablesSet::clone() const
{
    NotImplemented
    return nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const word& variablesSet::solverName() const
{
    return solverName_;
}


bool variablesSet::useSolverNameForFields() const
{
    return useSolverNameForFields_;
}


void variablesSet::setFluxField
(
    autoPtr<surfaceScalarField>& fieldPtr,
    const fvMesh& mesh,
    const volVectorField& velocity,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    // Try to read in field with custom or base name
    bool fieldFound
    (
        readFieldOK
        (
            fieldPtr,
            mesh,
            baseName,
            solverName,
            useSolverNameForFields
        )
    );

    // No base or custom field found.
    // Construct field based on linear interpolation
    if (!fieldFound)
    {
        word phiName(baseName);
        if (useSolverNameForFields)
        {
            phiName += solverName;
        }
        IOobject header
        (
            phiName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );
        fieldPtr.reset
        (
            new surfaceScalarField
            (
                header,
                linearInterpolate(velocity) & mesh.Sf()
            )
        );
    }
}


tmp<surfaceScalarField> variablesSet::allocateFluxField
(
    const fvMesh& mesh,
    const volVectorField& velocity,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    autoPtr<surfaceScalarField> fieldPtr(nullptr);
    setFluxField
    (
        fieldPtr,
        mesh,
        velocity,
        baseName,
        solverName,
        useSolverNameForFields
    );

    return tmp<surfaceScalarField>(fieldPtr.ptr());
}


tmp<volVectorField> variablesSet::autoCreateMeshMovementField
(
    const fvMesh& mesh,
    const word& fieldName,
    const dimensionSet& dims
)
{
    return tmp<volVectorField>::New
    (
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dims, Zero),
        fixedValueFvPatchVectorField::typeName
    );
}


void variablesSet::transfer(variablesSet& vars)
{
    // Does nothing in base
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
