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

#include "adjointSensitivity.H"
#include "sensitivityShapeFI.H"
#include "adjointSolver.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityShapeFI, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity, sensitivityShapeFI, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityShapeFI::sensitivityShapeFI
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    ShapeSensitivitiesBase(mesh, dict, adjointSolver)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivityShapeFI::computeDxDbInternalField() const
{
    return true;
}


void sensitivityShapeFI::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    if (eikonalSolver_)
    {
        eikonalSolver_->solve();
    }
    adjointSensitivity::assembleSensitivities(designVars);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
