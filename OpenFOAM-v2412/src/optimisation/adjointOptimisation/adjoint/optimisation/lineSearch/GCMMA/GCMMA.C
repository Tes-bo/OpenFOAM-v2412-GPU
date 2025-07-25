/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "GCMMA.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GCMMA, 1);
    addToRunTimeSelectionTable
    (
        lineSearch,
        GCMMA,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::GCMMA::writeToFiles(const bool isConverged)
{
    const PtrList<scalarField>& objValues = mma_.getValuesAndApproximations();
    const scalarField& rhoValues = mma_.getRho();
    const label m(rhoValues.size() - 1);
    if (Pstream::master())
    {
        unsigned int width = IOstream::defaultPrecision() + 5;
        if (writeHeader_)
        {
            GCMMAFile_
                << setw(width) << "#OuterIter" << " "
                << setw(width) << "InnerIter" << " "
                << setw(width) << "rhoObj" << " ";

            costFile_
                << setw(width) << "#nCycle" << " "
                << setw(width) << "cumulativeCost" << " "
                << setw(width) << "Objective" << " ";
            for (label i = 0; i < m; ++i)
            {
                GCMMAFile_
                    << setw(width) << "rhoConst" << " ";
                costFile_
                    << setw(width) << "Constraint" << " ";
            }
            GCMMAFile_
                << setw(width) << "J" << " "
                << setw(width) << "JTilda" << " ";
            for (label i = 0; i < m; ++i)
            {
                GCMMAFile_
                    << setw(width) << "C" << " "
                    << setw(width) << "CTilda" << " ";
            }
            GCMMAFile_<< endl;
            costFile_<< endl;

            writeHeader_ = false;
        }
        GCMMAFile_
            << setw(width) << iter_ + 2 << " "
            << setw(width) << innerIter_ << " ";

        forAll(rhoValues, i)
        {
            GCMMAFile_
                << setw(width) << rhoValues[i] << " ";
        }
        forAll(rhoValues, i)
        {
            GCMMAFile_
                << setw(width) << objValues[0][i] << " "
                << setw(width) << objValues[1][i] << " ";
        }
        GCMMAFile_ << endl;

        if (isConverged)
        {
            // The cost of this cycle is equal to the number of innerIters
            // plus one for the adjoint solution
            cost_ += innerIter_ + 1;
            costFile_
                << setw(width) << iter_ + 2 << " "
                << setw(width) << cost_ << " ";
            forAll(rhoValues, i)
            {
                costFile_
                    << setw(width) << objValues[0][i] << " ";
            }
            costFile_<< endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GCMMA::GCMMA
(
    const dictionary& dict,
    const Time& time,
    updateMethod& UpdateMethod
)
:
    lineSearch(dict, time, UpdateMethod),
    mma_(refCast<MMA>(UpdateMethod)),
    GCMMAFile_
        (time.globalPath()/"optimisation"/"objective"/time.timeName()/"GCMMA"),
    costFile_
    (
        time.globalPath()/"optimisation"/"objective"/time.timeName()
       /"GCMMACost"
    ),
    cost_(2),
    writeHeader_(true)
{
    // Force rho constants to be updated in each optimisation cycle
    mma_.setVariableRho();
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::GCMMA::converged()
{
    // Calls MMA::updateValuessAndApproximations() and needs to be called
    // before writeToFiles
    bool isConverged = mma_.converged();
    writeToFiles(isConverged);

    DebugInfo
        << "GCMMA converged ... " << Switch(isConverged) << endl;
    return isConverged;
}


void Foam::GCMMA::updateStep()
{
    DebugInfo
        << "GCMMA:: recomputing direction "<< endl;

    // Update rho values
    mma_.updateRho();

    // Re-solve the subproblem
    mma_.solveSubproblem();
}


void Foam::GCMMA::updateCorrection(scalarField& correction)
{
    correction = mma_.returnCorrection();
}


// ************************************************************************* //
