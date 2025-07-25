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

#include "SIMPLEControlOpt.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SIMPLEControlOpt, 1);
    addToRunTimeSelectionTable( SIMPLEControl, SIMPLEControlOpt, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::SIMPLEControlOpt::read()
{
    nInitialIters_ = dict().getOrDefault<label>("nInitialIters", nIters_);
    return SIMPLEControl::read();
}


bool Foam::SIMPLEControlOpt::criteriaSatisfied()
{
    bool satisfied(false);

    // Do not check criteria in the first iteration of the algorithm.
    // Used to avoid stopping the solution of the flow equations
    // due to a converged solution in the previous optimisation cycle
    if (subCycledTimePtr_().index() == 1)
    {
        satisfied = false;
    }
    else
    {
        satisfied = simpleControl::criteriaSatisfied();
    }

    return satisfied;
}


const Foam::label& Foam::SIMPLEControlOpt::nIters() const
{
    if (mesh_.time().timeIndex() == mesh_.time().startTimeIndex() + 1)
    {
        return nInitialIters_;
    }
    else
    {
        return nIters_;
    }
}


void Foam::SIMPLEControlOpt::resetDeltaT()
{
    Time& runTime = const_cast<Time&>(mesh_.time());
    if (runTime.deltaTValue() != deltaTSubSycle_)
    {
        runTime.setDeltaT(deltaTSubSycle_, false);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIMPLEControlOpt::SIMPLEControlOpt
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
:
    SIMPLEControl(mesh, managerType, solver),
    nInitialIters_(0),
    subCycledTimePtr_(nullptr),
    deltaTSubSycle_(Zero)
{
    this->read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SIMPLEControlOpt::write(const bool valid) const
{
    // Intentionally does nothing
    // Write is called only at the last iteration of the cycle
    return false;
}


bool Foam::SIMPLEControlOpt::loop()
{
    this->read();

    Time& runTime = const_cast<Time&>(mesh_.time());

    bool doNextIter(true);

    if (nIters() > 0)
    {
        // Sub-cycle time if this is the first iter
        if (!subCycledTimePtr_)
        {
            subCycledTimePtr_.reset(new subCycleTime(runTime, nIters()));
            Info<< "Solving equations for solver "
                << solver_.solverName() << "\n" << endl;
            deltaTSubSycle_ = runTime.deltaTValue();

            // Reset iteration count to zero
            iter_ = 0;

            // Reset previous time index of fvMesh::data, to avoid the rare
            // occurance of a solver satisfying the convergence criteria at the
            // first iteration, which then causes all subsequent optimisation
            // cycles to be seen as converged, irrespective of the residual
            // level, since the data::prevTimeIndex_ is not updated
            //mesh_.data::setPreviousTimeIndex(0);
        }

        // Increase index
        subCycledTimePtr_()++;
        iter_ = subCycledTimePtr_().index();


        if (criteriaSatisfied())
        {
            Info<< nl
                << solver_.solverName()
                << " solution converged in "
                << subCycledTimePtr_->index() << " iterations" << nl << endl;

            subCycledTimePtr_->endSubCycle();
            subCycledTimePtr_.clear();

            // Write solution before continuing to next solver
            runTime.write();
            solver_.write();

            // Check whether mean fields have not been computed due to an
            // unexpectedly early convergence
            checkMeanSolution();

            doNextIter = false;
        }
        else if (subCycledTimePtr_->end())
        {
            Info<< nl
                << solver_.solverName()
                << " solution reached max. number of iterations "
                << subCycledTimePtr_().nSubCycles() << nl << endl;

            subCycledTimePtr_->endSubCycle();
            subCycledTimePtr_.clear();

            // Write solution before continuing to next solver
            runTime.write();
            solver_.write();

            doNextIter = false;
        }
        else
        {
            // Since dicts are not updated when Time is sub-cycled,
            // do it manually here
            runTime.readModifiedObjects();
            resetDeltaT();

            DebugInfo
                << "Iteration " << subCycledTimePtr_().index()
                << "|" << subCycledTimePtr_().nSubCycles() << endl;

            simpleControl::storePrevIterFields();

            doNextIter = true;
        }
    }
    else
    {
        WarningInFunction
            << "Number of iterations is non-positive (" << nIters() << ")."
            << nl
            << "Skipping the solution of the equations corresponding to solver "
            << solver_.solverName()
            << nl << endl;

        doNextIter = false;
    }

    return doNextIter;
}


// ************************************************************************* //
