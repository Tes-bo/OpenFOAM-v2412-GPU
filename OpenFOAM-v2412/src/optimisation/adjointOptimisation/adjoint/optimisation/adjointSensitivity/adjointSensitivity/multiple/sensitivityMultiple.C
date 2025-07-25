/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "adjointSensitivity.H"
#include "sensitivityMultiple.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityMultiple, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivityMultiple,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityMultiple::sensitivityMultiple
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    adjointSensitivity(mesh, dict, adjointSolver),
    sensTypes_(this->dict().get<wordList>("sensitivityTypes")),
    sens_(sensTypes_.size())
{
    forAll(sensTypes_, sI)
    {
        sens_.set
        (
            sI,
            adjointSensitivity::New
            (
                mesh,
                this->dict().subDict(sensTypes_[sI]),
                adjointSolver
            )
        );
        sens_[sI].setSuffix(sensTypes_[sI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivityMultiple::readDict(const dictionary& dict)
{
    if (adjointSensitivity::readDict(dict))
    {
        forAll(sens_, sI)
        {
            sens_[sI].readDict(dict.subDict(sensTypes_[sI]));
        }

        return true;
    }

    return false;
}


void sensitivityMultiple::accumulateIntegrand(const scalar dt)
{
    forAll(sens_, sI)
    {
        sens_[sI].accumulateIntegrand(dt);
    }
}


void sensitivityMultiple::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    forAll(sens_, sI)
    {
        sens_[sI].assembleSensitivities(designVars);
    }
}


const scalarField& sensitivityMultiple::calculateSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    forAll(sens_, sI)
    {
        Info<< "Computing sensitivities " << sensTypes_[sI] << endl;
        sens_[sI].calculateSensitivities(designVars);
    }
    write(type());

    return derivatives_;
}


void sensitivityMultiple::clearSensitivities()
{
    forAll(sens_, sI)
    {
        sens_[sI].clearSensitivities();
    }
}


void sensitivityMultiple::write(const word& baseName)
{
    forAll(sens_, sI)
    {
        sens_[sI].write(sensTypes_[sI]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
