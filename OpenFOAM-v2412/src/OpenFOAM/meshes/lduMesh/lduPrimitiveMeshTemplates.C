/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "lduPrimitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ProcPatch>
Foam::lduSchedule Foam::lduPrimitiveMesh::nonBlockingSchedule
(
    const lduInterfacePtrsList& interfaces
)
{
    lduSchedule schedule(2*interfaces.size());

    label patchEvali = 0;
    label numProcPatches = 0;

    //
    // 1. Schedule non-processor patches
    //

    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            if (isA<ProcPatch>(interfaces[patchi]))
            {
                ++numProcPatches;
            }
            else
            {
                schedule[patchEvali++].setInitEvaluate(patchi);
                schedule[patchEvali++].setEvaluate(patchi);
            }
        }
    }


    //
    // 2. Schedule processor patches
    //

    if (numProcPatches)
    {
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi) && isA<ProcPatch>(interfaces[patchi]))
            {
                schedule[patchEvali].setInitEvaluate(patchi);
                schedule[patchEvali + numProcPatches].setEvaluate(patchi);
                ++patchEvali;
            }
        }
    }

    // Caution:
    // The schedule is only valid for a subset of its range
    // (where interfaces are defined) but must retain the full list length
    // for later (external) bookkeeping

    schedule.setSize(patchEvali);

    return schedule;
}


// ************************************************************************* //
