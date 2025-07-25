/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "RASTurbulenceModel.H"
#include "findRefCell.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RASTurbulenceModel, 0);
    addToRunTimeSelectionTable
    (
        incompressiblePrimalSolver,
        RASTurbulenceModel,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::incompressibleVars& Foam::RASTurbulenceModel::allocateVars()
{
    vars_.reset(new incompressibleVars(mesh_, solverControl_()));
    return getIncoVars();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASTurbulenceModel::RASTurbulenceModel
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& solverName
)
:
    incompressiblePrimalSolver
    (
        mesh,
        managerType,
        dict,
        solverName
    ),
    solverControl_(SIMPLEControl::New(mesh, managerType, *this)),
    incoVars_(allocateVars())
{
    setRefCell
    (
        incoVars_.pInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RASTurbulenceModel::solveIter()
{
    const Time& time = mesh_.time();
    Info<< "Time = " << time.timeName() << "\n" << endl;

    // Grab references
    autoPtr<incompressible::turbulenceModel>& turbulence =
        incoVars_.turbulence();
    turbulence->correct();

    solverControl_().write();

    // Average fields if necessary
    incoVars_.computeMeanFields();

    time.printExecutionTime(Info);
}


void Foam::RASTurbulenceModel::solve()
{
    // Iterate
    if (active_)
    {
        // Reset initial and mean fields before solving
        while (solverControl_().loop())
        {
            solveIter();
        }
    }
}


bool Foam::RASTurbulenceModel::loop()
{
    return solverControl_().loop();
}


bool Foam::RASTurbulenceModel::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}


// ************************************************************************* //
