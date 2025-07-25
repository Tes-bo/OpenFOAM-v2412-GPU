/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "objectiveFlowRate.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveFlowRate, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveFlowRate,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveFlowRate::objectiveFlowRate
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    patches_
    (
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches")).sortedToc()
    ),
    flowRates_(patches_.size(), Zero)
{
    // Allocate boundary field pointers
    bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveFlowRate::J()
{
    J_ = 0;
    const surfaceScalarField& phi = vars_.phiInst();

    forAll(patches_, pI)
    {
        const label patchI = patches_[pI];
        flowRates_[pI] = gSum(phi.boundaryField()[patchI]);
        J_ += flowRates_[pI];
    }

    return J_;
}


void objectiveFlowRate::update_boundarydJdv()
{
    for (const label patchI : patches_)
    {
        bdJdvPtr_()[patchI] = mesh_.boundary()[patchI].nf();
    }
}


void objectiveFlowRate::update_boundarydJdvn()
{
    for (const label patchI : patches_)
    {
        bdJdvnPtr_()[patchI] = 1;
    }
}


void objectiveFlowRate::addHeaderColumns() const
{
    for (const label patchI : patches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        objFunctionFilePtr_()
            << setw(width_) << word(patch.name() + "FlowRate") << " ";
    }
}


void objectiveFlowRate::addColumnValues() const
{
    for (const scalar flowRate : flowRates_)
    {
        objFunctionFilePtr_()
            << setw(width_) << flowRate << " ";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
