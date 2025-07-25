/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "SIMPLEControlSingleRun.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SIMPLEControlSingleRun, 0);
    addToRunTimeSelectionTable
    (
        SIMPLEControl,
        SIMPLEControlSingleRun,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::SIMPLEControlSingleRun::read()
{
    return SIMPLEControl::read();
}


void Foam::SIMPLEControlSingleRun::readIters()
{
    label nItersOld = nIters_;
    nIters_ = dict().get<label>("nIters");

    if (nIters_ != nItersOld || iter_ == 0)
    {
        Time& runTime = const_cast<Time&>(mesh_.time());
        if (iter_ == 0)
        {
            startTime_ = runTime.value();
        }
        Info<< "Setting endTime to " << startTime_ + nIters_ << endl;
        runTime.setEndTime(startTime_ + nIters_);
        endTime_ = runTime.endTime().value();
    }
}


void Foam::SIMPLEControlSingleRun::checkEndTime(bool& isRunning)
{
    // If controlDict is modified during run-time, time.endTime() is reset
    // to what is read from controlDict and overwrites the one set through
    // nIters. Silently reset
    Time& time = const_cast<Time&>(mesh_.time());

    if (time.endTime().value() != endTime_)
    {
        time.setEndTime(startTime_ + nIters_);
        endTime_ = time.endTime().value();
        isRunning =
            time.value() < (time.endTime().value() - 0.5*time.deltaTValue());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIMPLEControlSingleRun::SIMPLEControlSingleRun
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
:
    SIMPLEControl(mesh, managerType, solver),
    startTime_(Zero),
    endTime_(Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SIMPLEControlSingleRun::write(const bool valid) const
{
    Time& time = const_cast<Time&>(mesh_.time());
    time.write();
    solver_.write();

    return true;
}


void Foam::SIMPLEControlSingleRun::writeNow()
{
    Time& time = const_cast<Time&>(mesh_.time());
    // Avoid writing fields if already in an writeTime iter
    // since results will be written by the solver class either way
    if (!time.writeTime())
    {
        time.writeNow();
        solver_.writeNow();
    }
}


bool Foam::SIMPLEControlSingleRun::loop()
{
    solutionControl::setFirstIterFlag(true, true);

    this->read();
    ++iter_;

    Time& runTime = const_cast<Time&>(mesh_.time());

    if (initialised_ && criteriaSatisfied())
    {
        Info<< nl
            << solver_.solverName()
            << " solution converged in "
            << runTime.timeName() << " iterations" << nl << endl;

        // write fields (including dummy turbulence fields in multi-point runs)
        writeNow();

        // Check whether mean fields have not been computed due to an
        // unexpectedly early convergence
        checkMeanSolution();

        return false;
    }
    else
    {
        initialised_ = true;
        simpleControl::storePrevIterFields();
    }

    bool isRunning = runTime.loop();
    checkEndTime(isRunning);

    if (!isRunning)
    {
        Info<< nl
            << solver_.solverName()
            << " solution reached max. number of iterations "
            << nIters_ << nl << endl;

        // Write fields (including dummy turbulence fields in multi-point runs)
        writeNow();
    }

    return isRunning;
}


// ************************************************************************* //
