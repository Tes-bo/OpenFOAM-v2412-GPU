/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::pimpleControl::read()
{
    solutionControl::read(false);

    const dictionary pimpleDict(dict());

    solveFlow_ = pimpleDict.getOrDefault("solveFlow", true);
    nCorrPIMPLE_ = pimpleDict.getOrDefault<label>("nOuterCorrectors", 1);
    nCorrPISO_ = pimpleDict.getOrDefault<label>("nCorrectors", 1);
    SIMPLErho_ = pimpleDict.getOrDefault("SIMPLErho", false);
    turbOnFinalIterOnly_ =
        pimpleDict.getOrDefault("turbOnFinalIterOnly", true);
    finalOnLastPimpleIterOnly_ =
        pimpleDict.getOrDefault("finalOnLastPimpleIterOnly", false);
    ddtCorr_ = pimpleDict.getOrDefault("ddtCorr", true);

    return true;
}


bool Foam::pimpleControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    const bool storeIni = this->storeInitialResiduals();

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

            scalar relative = 0.0;
            bool relCheck = false;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = residuals.first();
            }
            else
            {
                const scalar iniRes =
                    (residualControl_[fieldi].initialResidual + ROOTVSMALL);

                relative = residuals.last() / iniRes;
                relCheck = (relative < residualControl_[fieldi].relTol);
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << fieldName
                    << " PIMPLE iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldi].initialResidual
                    << ", abs tol = " << residuals.last()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldi].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::pimpleControl::setFirstIterFlag(const bool check, const bool force)
{
    DebugInfo
        << "corr:" << corr_
        << " corrPISO:" << corrPISO_
        << " corrNonOrtho:" << corrNonOrtho_
        << endl;

    solutionControl::setFirstIterFlag(check && corrPISO_ <= 1, force);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleControl::pimpleControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool verbose
)
:
    solutionControl(mesh, dictName),
    solveFlow_(true),
    nCorrPIMPLE_(0),
    nCorrPISO_(0),
    corrPISO_(0),
    SIMPLErho_(false),
    turbOnFinalIterOnly_(true),
    finalOnLastPimpleIterOnly_(false),
    ddtCorr_(true),
    converged_(false)
{
    read();

    if (verbose)
    {
        Info<< nl << algorithmName_;

        if (nCorrPIMPLE_ > 1)
        {
            if (residualControl_.empty())
            {
                Info<< ": no residual control data found. "
                    << "Calculations will employ " << nCorrPIMPLE_
                    << " corrector loops" << nl;
            }
            else
            {
                Info<< ": max iterations = " << nCorrPIMPLE_ << nl;

                for (const fieldData& ctrl : residualControl_)
                {
                    Info<< "    field " << ctrl.name << token::TAB
                        << ": relTol " << ctrl.relTol
                        << ", tolerance " << ctrl.absTol
                        << nl;
                }
            }
        }
        else
        {
            Info<< ": Operating solver in PISO mode" << nl;
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pimpleControl::loop()
{
    read();

    ++corr_;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    setFirstIterFlag();

    if (corr_ == nCorrPIMPLE_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPIMPLE_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data().setFinalIteration(false);
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;

            mesh_.data().setFinalIteration(false);
            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();

            mesh_.data().setFinalIteration(true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data().setFinalIteration(true);
        }

        if (corr_ <= nCorrPIMPLE_)
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
