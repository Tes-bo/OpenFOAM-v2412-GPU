/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "adjointSensitivity.H"
#include "sensitivityTopO.H"
#include "adjointSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityTopO, 0);
addToRunTimeSelectionTable(adjointSensitivity, sensitivityTopO, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void sensitivityTopO::zeroSensInFixedPorousZones(scalarField& sens)
{
    const labelList& adjointPorousIDs = zones_.adjointPorousZoneIDs();
    if (adjointPorousIDs.empty())
    {
        for (label cellZoneID : zones_.fixedPorousZoneIDs())
        {
            const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
            for (label cellI : zoneCells)
            {
                sens[cellI] = 0.;
            }
        }
        for (label cellZoneID : zones_.fixedZeroPorousZoneIDs())
        {
            const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
            for (label cellI : zoneCells)
            {
                sens[cellI] = 0.;
            }
        }
        for (label cellI : zones_.IOCells())
        {
            sens[cellI] = 0.;
        }
    }
    else
    {
        // if adjointPorousZones are not empty, zero sensitivities in all cells
        // not within them
        scalarField mask(sens.size(), Zero);
        for (label cellZoneID : adjointPorousIDs)
        {
            const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
            for (label cellI : zoneCells)
            {
                mask[cellI] = 1.;
            }
        }
        sens *= mask;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityTopO::sensitivityTopO
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    adjointSensitivity(mesh, dict, adjointSolver),
    zones_(mesh, dict.parent()),
    designVariablesName_("beta")
{
    if (includeDistance_)
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict_,
                adjointSolver,
                labelHashSet(0)
            )
        );
    }
    // Allocate sensitivities field
    fieldSensPtr_.reset
    (
        createZeroFieldPtr<scalar>
        (
            mesh_,
            "topologySens" + adjointSolver.solverName(),
            pow5(dimLength)/sqr(dimTime)
        )
    );

    // Set return field size
    derivatives_ = scalarField(mesh_.nCells(), Zero);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivityTopO::readDict(const dictionary& dict)
{
    if (adjointSensitivity::readDict(dict))
    {
        if (includeDistance_)
        {
            if (eikonalSolver_)
            {
                eikonalSolver_->readDict(dict);
            }
            else
            {
                eikonalSolver_.reset
                (
                    new adjointEikonalSolver
                    (
                        mesh_,
                        dict_,
                        adjointSolver_,
                        labelHashSet(0)
                    )
                );
            }
        }

        return true;
    }

    return false;
}


void sensitivityTopO::accumulateIntegrand(const scalar dt)
{
    // Accumulate source for additional post-processing PDEs, if necessary
    if (eikonalSolver_)
    {
        eikonalSolver_->accumulateIntegrand(dt);
    }

    adjointSolver_.topOSensMultiplier
        (fieldSensPtr_().primitiveFieldRef(), designVariablesName_, dt);
}


void sensitivityTopO::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    scalarField& sens = fieldSensPtr_().primitiveFieldRef();
    if (eikonalSolver_)
    {
        eikonalSolver_->solve();
        sens += eikonalSolver_->topologySensitivities(designVariablesName_);
    }
    zeroSensInFixedPorousZones(sens);

    adjointSensitivity::assembleSensitivities(designVars);
}


void sensitivityTopO::postProcessSens
(
    scalarField& sens,
    scalarField& auxSens,
    fv::options& fvOptions,
    const word& fieldName,
    const word& designVariablesName
)
{
    if (fvOptions.appliesToField(fieldName))
    {
        DebugInfo
            << "Computing SD contributions from the interpolation of "
            << fieldName << endl;
        fvOptions.postProcessSens(auxSens, fieldName, designVariablesName);
        sens += auxSens;
    }
}


void sensitivityTopO::postProcessSens
(
    scalarField& sens,
    scalarField& auxSens,
    const word& fieldName
)
{
    fv::options& fvOptions(fv::options::New(this->mesh_));
    postProcessSens(sens, auxSens, fvOptions, fieldName, designVariablesName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
