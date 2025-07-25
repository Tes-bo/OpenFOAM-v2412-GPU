/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "simpleControl.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::simpleControl::read()
{
    solutionControl::read(true);
    return true;
}


bool Foam::simpleControl::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.data().solverPerformanceDict();
    for (const entry& solverPerfDictEntry : solverDict)
    {
        const word& fieldName = solverPerfDictEntry.keyword();
        const label fieldi = applyToField(fieldName);

        if (fieldi != -1)
        {
            Pair<scalar> residuals = maxResidual(solverPerfDictEntry);
            checked = true;

            const bool absCheck =
                (residuals.first() < residualControl_[fieldi].absTol);

            achieved = achieved && absCheck;

            if (debug)
            {
                Info<< algorithmName_ << " solution statistics:" << endl;

                Info<< "    " << fieldName << ": tolerance = "
                    << residuals.first()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleControl::simpleControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool verbose
)
:
    solutionControl(mesh, dictName),
    initialised_(false)
{
    read();

    if (verbose)
    {
        Info<< nl << algorithmName_;

        if (residualControl_.empty())
        {
            const scalar duration =
                mesh_.time().endTime().value()
              - mesh_.time().startTime().value();

            Info<< ": no convergence criteria found. "
                << "Calculations will run for " << duration << " steps."
                << nl;
        }
        else
        {
            Info<< ": convergence criteria" << nl;
            for (const fieldData& ctrl : residualControl_)
            {
                Info<< "    field " << ctrl.name << token::TAB
                    << " tolerance " << ctrl.absTol
                    << nl;
            }
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simpleControl::loop()
{
    solutionControl::setFirstIterFlag(true, true);

    read();

    Time& runTime = const_cast<Time&>(mesh_.time());

    if (initialised_ && criteriaSatisfied())
    {
        Info<< nl
            << algorithmName_
            << " solution converged in "
            << runTime.timeName() << " iterations" << nl << endl;

        // Set to finalise calculation
        runTime.writeAndEnd();
    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }

    return runTime.loop();
}


// ************************************************************************* //
