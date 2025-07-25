/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "orderedPhasePair.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orderedPhasePair::orderedPhasePair
(
    const multiphaseInter::phaseModel& from,
    const multiphaseInter::phaseModel& to
)
:
    phasePair(from, to, true)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::multiphaseInter::phaseModel& Foam::orderedPhasePair::from() const
{
    return phase1();
}


const Foam::multiphaseInter::phaseModel& Foam::orderedPhasePair::to() const
{
    return phase2();
}


Foam::word Foam::orderedPhasePair::name() const
{
    word namec(second());
    namec[0] = toupper(namec[0]);
    return first() + "To" + namec;
}


// ************************************************************************* //
