/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::lduScheduleEntry

Description
    Struct to hold the patch index and the initialisation flag for
    the patch schedule

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduSchedule_H
#define Foam_lduSchedule_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

struct lduScheduleEntry;
Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb);

//- A List of lduSchedule entries.
typedef List<lduScheduleEntry> lduSchedule;


/*---------------------------------------------------------------------------*\
                      Class lduScheduleEntry Declaration
\*---------------------------------------------------------------------------*/

struct lduScheduleEntry
{
    label patch;
    bool init;

    //- Define patch for initEvaluate()
    void setInitEvaluate(const label patchi) noexcept
    {
        patch = patchi;
        init = true;
    }

    //- Define patch for evaluate()
    void setEvaluate(const label patchi) noexcept
    {
        patch = patchi;
        init = false;
    }

    //- Always unequal: each entry is supposed to be unique!
    friend bool operator!=(const lduScheduleEntry&, const lduScheduleEntry&)
    {
        return true;
    }

    friend Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb)
    {
        os  << lb.patch << token::SPACE << lb.init;
        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
