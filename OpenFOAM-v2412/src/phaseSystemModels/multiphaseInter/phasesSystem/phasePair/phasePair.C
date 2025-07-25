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

#include "phasePair.H"
#include "surfaceTensionModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePair::phasePair
(
    const multiphaseInter::phaseModel& phase1,
    const multiphaseInter::phaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::multiphaseInter::phaseModel& Foam::phasePair::from() const
{
    FatalErrorInFunction
        << "Requested from phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::multiphaseInter::phaseModel& Foam::phasePair::to() const
{
    FatalErrorInFunction
        << "Requested to phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


Foam::word Foam::phasePair::name() const
{
    word name2(second());
    name2[0] = toupper(name2[0]);
    return first() + "And" + name2;
}


// ************************************************************************* //
