/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "fieldExprDriver.H"
#include "fieldExprScanner.H"
#include "error.H"
#include "className.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace fieldExpr
{

defineTypeNameAndDebug(parseDriver, 0);

} // End namespace fieldExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::fieldExpr::parseDriver::parseDriver
(
    const label len
)
:
    parsing::genericRagelLemonDriver(),
    expressions::exprDriver(),
    size_(Foam::max(1, len))
{}


Foam::expressions::fieldExpr::parseDriver::parseDriver
(
    const label len,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::exprDriver(dict),
    size_(len)
{}


Foam::expressions::fieldExpr::parseDriver::parseDriver
(
    const label len,
    const parseDriver& rhs,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::exprDriver(rhs, dict),
    size_(len)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

unsigned Foam::expressions::fieldExpr::parseDriver::parse
(
    const std::string& expr,
    size_t pos,
    size_t len
)
{
    scanner scan(this->debugScanner());

    scan.process(expr, pos, len, *this);

    return 0;
}


// ************************************************************************* //
