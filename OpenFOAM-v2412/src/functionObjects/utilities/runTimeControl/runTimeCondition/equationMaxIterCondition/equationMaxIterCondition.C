/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "equationMaxIterCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(equationMaxIterCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        equationMaxIterCondition,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::equationMaxIterCondition::
equationMaxIterCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    fieldNames_(dict.get<wordList>("fields")),
    threshold_(dict.get<label>("threshold")),
    startIter_(dict.getOrDefault("startIter", 2))
{
    if (!fieldNames_.size())
    {
        WarningInFunction
            << "No fields supplied: deactivating" << endl;

        active_ = false;
    }

    startIter_ = max(startIter_, 2);
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::equationMaxIterCondition::apply()
{
    bool satisfied = false;

    if (!active_)
    {
        return true;
    }

    if (obr_.time().timeIndex() < startIter_)
    {
        // Do not start checking until start iter
        return false;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const dictionary& solverDict = mesh.data().solverPerformanceDict();

    List<label> result(fieldNames_.size(), -1);

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        if (solverDict.found(fieldName))
        {
            const List<solverPerformance> sp(solverDict.lookup(fieldName));
            const label nIterations = sp.first().nIterations();
            result[fieldi] = nIterations;

            if (nIterations > threshold_)
            {
                satisfied = true;
            }
        }
    }

    bool valid = false;
    forAll(result, i)
    {
        if (result[i] < 0)
        {
            WarningInFunction
                << "Number of iterations data not found for field "
                << fieldNames_[i] << endl;
        }
        else
        {
            valid = true;
        }
    }

    if (!valid)
    {
        WarningInFunction
            << "Number of iterations data not found for any fields: "
            << "deactivating" << endl;

        active_ = false;
    }

    if (satisfied && valid)
    {
        Log << type() << ": " << name_
            << ": satisfied using threshold value: " << threshold_ << nl;

        forAll(result, resulti)
        {
            if (result[resulti] != -1)
            {
                Log << "    field: " << fieldNames_[resulti]
                    << ", iterations: " << result[resulti] << nl;
            }
        }
        Log << endl;
    }

    return satisfied;
}


void Foam::functionObjects::runTimeControls::equationMaxIterCondition::write()
{
    // do nothing
}



void Foam::functionObjects::runTimeControls::equationMaxIterCondition::reset()
{
    // do nothing
}


// ************************************************************************* //
