/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
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

#include "exprResultStoredStack.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

    defineTypeName(exprResultStoredStack);

    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStoredStack,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStoredStack,
        empty
    );

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprResultStoredStack::exprResultStoredStack()
:
    expressions::exprResultStack()
{
    needsReset(false);  // Override parent: does not reset every timestep
}


Foam::expressions::exprResultStoredStack::exprResultStoredStack
(
    const exprResultStoredStack& rhs
)
:
    expressions::exprResultStack(rhs)
{
    needsReset(false);  // Override parent: does not reset every timestep
}


Foam::expressions::exprResultStoredStack::exprResultStoredStack
(
    const dictionary& dict
)
:
    expressions::exprResultStack(dict)
{
    needsReset(false);  // Override parent: does not reset every timestep
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::expressions::exprResultStoredStack::operator=
(
    const exprResultStoredStack& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    exprResultStack::operator=(rhs);
}


// ************************************************************************* //
