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

#include "objectiveFlowRatePartition.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveFlowRatePartition, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveFlowRatePartition,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveFlowRatePartition::objectiveFlowRatePartition
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    inletPatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("inletPatches")
        ).sortedToc()
    ),
    outletPatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("outletPatches")
        ).sortedToc()
    ),
    targetFlowRateFraction_(),
    currentFlowRateFraction_(outletPatches_.size(), Zero),
    inletFlowRate_(0),
    flowRateDifference_(outletPatches_.size(), Zero)
{
    // Read target fractions if present, otherwise treat them as equally
    // distributed
    if
    (
        !dict.readIfPresentCompat
        (
            "targetFractions", {{"targetPercentages", 2306}},
            targetFlowRateFraction_
        )
    )
    {
        const label nOutPatches = outletPatches_.size();
        targetFlowRateFraction_.setSize(nOutPatches, 1.0/scalar(nOutPatches));
    }
    // Sanity checks
    if (targetFlowRateFraction_.size() != outletPatches_.size())
    {
        FatalErrorInFunction
            << "Inconsistent sizes for targetFractions and outletPatches"
            << exit(FatalError);

    }

    // Allocate boundary field pointers
    bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));

    // Determine word width for output
    for (const label patchI : outletPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        unsigned int wordSize = word(patch.name() + "Ratio").size();
        width_ = max(width_, wordSize);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveFlowRatePartition::J()
{
    J_ = 0;
    const surfaceScalarField& phi = vars_.phiInst();

    // Inlet patches
    inletFlowRate_ = 0;
    for (const label patchI : inletPatches_)
    {
        // Negative value
        inletFlowRate_ += gSum(phi.boundaryField()[patchI]);
    }

    // Outlet patches
    forAll(outletPatches_, pI)
    {
        const label patchI = outletPatches_[pI];
        const scalar outletFlowRate = gSum(phi.boundaryField()[patchI]);
        currentFlowRateFraction_[pI] = -outletFlowRate/inletFlowRate_;
        flowRateDifference_[pI] =
            targetFlowRateFraction_[pI] - currentFlowRateFraction_[pI];
        J_ += 0.5*flowRateDifference_[pI]*flowRateDifference_[pI];
    }

    return J_;
}


void objectiveFlowRatePartition::update_boundarydJdv()
{
    forAll(outletPatches_, pI)
    {
        const label patchI = outletPatches_[pI];
        tmp<vectorField> nf = mesh_.boundary()[patchI].nf();
        bdJdvPtr_()[patchI] = nf*flowRateDifference_[pI]/inletFlowRate_;
    }
}


void objectiveFlowRatePartition::update_boundarydJdvn()
{
    forAll(outletPatches_, pI)
    {
        const label patchI = outletPatches_[pI];
        bdJdvnPtr_()[patchI] = flowRateDifference_[pI]/inletFlowRate_;
    }
}


void objectiveFlowRatePartition::addHeaderInfo() const
{
    objFunctionFilePtr_()
        << setw(width_) << "#inletFlowRate" << " "
        << setw(width_) << -inletFlowRate_ << endl;
    forAll(outletPatches_, pI)
    {
        const label patchI = outletPatches_[pI];
        const fvPatch& patch = mesh_.boundary()[patchI];
        objFunctionFilePtr_()
            << setw(width_) << word("#" + patch.name() + "Tar") << " "
            << setw(width_) << targetFlowRateFraction_[pI] << endl;
    }
}


void objectiveFlowRatePartition::addHeaderColumns() const
{
    for (const label patchI : outletPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        objFunctionFilePtr_()
            << setw(width_) << word(patch.name() + "Ratio") << " ";
    }
}


void objectiveFlowRatePartition::addColumnValues() const
{
    for (const scalar flowRate : currentFlowRateFraction_)
    {
        objFunctionFilePtr_()
            << setw(width_) << flowRate << " ";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
