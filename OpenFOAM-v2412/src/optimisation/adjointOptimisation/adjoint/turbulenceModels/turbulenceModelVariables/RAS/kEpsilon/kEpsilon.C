/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "kEpsilon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilon, 0);
addToRunTimeSelectionTable(RASModelVariables, kEpsilon, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kEpsilon::allocateMeanFields()
{
    RASModelVariables::allocateMeanFields();
    if (solverControl_.average())
    {
        GMean_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "GMean",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimArea/pow3(dimTime), Zero)
            )
        );
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilon::kEpsilon
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    RASModelVariables(mesh, SolverControl)
{
    TMVar1BaseName_ = "k";
    TMVar2BaseName_ = "epsilon";

    TMVar1Ptr_.ref(mesh_.lookupObjectRef<volScalarField>(TMVar1BaseName_));
    TMVar2Ptr_.ref(mesh_.lookupObjectRef<volScalarField>(TMVar2BaseName_));
    nutPtr_.ref(mesh_.lookupObjectRef<volScalarField>(nutBaseName_));

    allocateInitValues();
    allocateMeanFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField::Internal> kEpsilon::computeG()
{
    const turbulenceModel& turbModel = mesh_.lookupObject<turbulenceModel>
    (
         IOobject::groupName
         (
             turbulenceModel::propertiesName,
             TMVar2().internalField().group()
         )
    );
    // Recompute G and modify values next to the walls
    // Ideally, grad(U) should be cached to avoid the overhead
    const volVectorField& U = turbModel.U();
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal GbyNu0
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        (tgradU() && devTwoSymm(tgradU()))
    );

    // NB: leave tmp registered (for correctBoundaryConditions)
    auto tG =
        tmp<volScalarField::Internal>::New
        (
            turbModel.GName(),
            nutRefInst()*GbyNu0
        );
    // Use correctBoundaryConditions instead of updateCoeffs to avoid
    // messing with updateCoeffs in the next iteration of omegaEqn
    TMVar2Inst().correctBoundaryConditions();

    return tG;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField::Internal> kEpsilon::G()
{
    if (solverControl_.useAveragedFields())
    {
        DebugInfo
            << "Using GMean" << endl;
        return tmp<volScalarField::Internal>(GMean_());
    }
    DebugInfo
        << "Using instantaneous G" << endl;
    return computeG();
}


void kEpsilon::computeMeanFields()
{
    RASModelVariables::computeMeanFields();
    if (solverControl_.doAverageIter())
    {
        const label iAverageIter = solverControl_.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter + 1);
        scalar mult = avIter*oneOverItP1;
        GMean_() = GMean_()*mult + computeG()*oneOverItP1;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
