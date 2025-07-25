/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "readFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::label Foam::checkData
(
    const objectRegistry& obr,
    const instantList& timeDirs,
    wordList& objectNames,
    const fileName& local
)
{
    // Assume prune_0() was used prior to calling this

    wordHashSet goodFields;

    for (const word& fieldName : objectNames)
    {
        // // If prune_0() not previously used...
        // if (objectNames.ends_with("_0")) continue;

        bool good = false;

        for (const instant& inst : timeDirs)
        {
            good =
                IOobject
                (
                    fieldName,
                    inst.name(),
                    local,
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ).typeHeaderOk<volScalarField>(false, false);

            if (!good)
            {
                break;
            }
        }

        if (returnReduceAnd(good))
        {
            goodFields.insert(fieldName);
        }
    }

    objectNames = goodFields.sortedToc();

    return objectNames.size();
}


// ************************************************************************* //
