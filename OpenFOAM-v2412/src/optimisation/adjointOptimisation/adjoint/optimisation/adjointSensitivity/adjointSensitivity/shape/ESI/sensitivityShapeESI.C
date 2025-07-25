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
    ESITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundaryFieldsFwd.H"
#include "sensitivityShapeESI.H"
#include "adjointSolver.H"
#include "ShapeSensitivitiesBase.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityShapeESI, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity, sensitivityShapeESI, dictionary
);


void sensitivityShapeESI::computeDxDbMult()
{
    if (eikonalSolver_)
    {
        eikonalSolver_->solve();
    }
    if (adjointMeshMovementSolver_)
    {
        adjointMeshMovementSolver_->solve();
        boundaryVectorField& meshMovementSens =
            adjointMeshMovementSolver_->meshMovementSensitivities();
        PtrList<objective>& functions =
            adjointSolver_.getObjectiveManager().getObjectiveFunctions();
        for (const label patchI : geometryVariationIntegrationPatches())
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            const scalarField& magSf = patch.magSf();
            const vectorField& Sf = patch.Sf();
            dxdbMult_()[patchI] = meshMovementSens[patchI]*magSf;
            for (objective& func : functions)
            {
                if (func.hasDivDxDbMult())
                {
                    Info<< func.objectiveName() << " " << patch.name() << endl;
                    dxdbDirectMult_()[patchI] +=
                        func.weight()
                       *func.divDxDbMultiplier().boundaryField()[patchI]
                       *Sf;
                }
            }
        }
    }
    for (const label patchI : geometryVariationIntegrationPatches())
    {
        const vectorField& Sf = mesh_.boundary()[patchI].Sf();
        dxdbMult_()[patchI] += Sf & gradDxDbMult_().boundaryField()[patchI];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityShapeESI::sensitivityShapeESI
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    ShapeSensitivitiesBase(mesh, dict, adjointSolver)
{
    dxdbMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    // The boundary values of divDxDbMultiplier are stored in dxdbDirectMult
    // after applying the Gauss divergence theorem.
    // Allocate dxdbDirectMult if necessary
    if (hasMultiplier(&objective::hasDivDxDbMult))
    {
        dxdbDirectMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    if (dict.getOrDefault<bool>("includeMeshMovement", true))
    {
        adjointMeshMovementSolver_.reset
        (
            new adjointMeshMovementSolver(mesh, dict, *this)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sensitivityShapeESI::readDict(const dictionary& dict)
{
    if (ShapeSensitivitiesBase::readDict(dict))
    {
        bool includeMeshMovement =
            dict.getOrDefault<bool>("includeMeshMovement", true);

        if (includeMeshMovement)
        {
            if (adjointMeshMovementSolver_)
            {
                adjointMeshMovementSolver_->readDict(dict);
            }
            else
            {
                adjointMeshMovementSolver_.reset
                (
                    new adjointMeshMovementSolver(mesh_, dict, *this)
                );
            }
        }

        return true;
    }

    return false;
}


void sensitivityShapeESI::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    computeDxDbMult();
    if (designVars)
    {
        adjointSensitivity::assembleSensitivities(designVars);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
