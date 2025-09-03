/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "exprDriver.H"
#include "FieldOps.H"
#include "Random.H"
#include "TimeState.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprDriver::fill_random
(
    scalarField& field,
    label seed,
    const bool gaussian
) const
{
    if (seed <= 0)
    {
        const TimeState* ts = this->timeState();

        if (ts)
        {
            seed = (ts->timeIndex() - seed);
        }
        else
        {
            seed = -seed;
        }
    }

    if (gaussian)
    {
        Random::gaussianGeneratorOp<scalar> gen(seed);
        FieldOps::assign(field, field, gen);
    }
    else
    {
        Random::uniformGeneratorOp<scalar> gen(seed);
        FieldOps::assign(field, field, gen);
    }
}


Foam::point Foam::expressions::exprDriver::getPositionOfMinimum
(
    const scalarField& vals,
    const pointField& locs
)
{
    return FieldOps::findMinData(vals, locs).second();
}


Foam::point Foam::expressions::exprDriver::getPositionOfMaximum
(
    const scalarField& vals,
    const pointField& locs
)
{
    return FieldOps::findMaxData(vals, locs).second();
}


// ************************************************************************* //
