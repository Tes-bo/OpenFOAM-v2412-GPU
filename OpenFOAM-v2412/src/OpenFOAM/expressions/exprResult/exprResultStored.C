/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "exprResultStored.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
    defineTypeName(exprResultStored);

    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStored,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStored,
        empty
    );

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprResultStored::exprResultStored()
:
    expressions::exprResult(),
    name_("none"),
    startExpr_()
{}


Foam::expressions::exprResultStored::exprResultStored
(
    const exprResultStored& rhs
)
:
    expressions::exprResult(rhs),
    name_(rhs.name_),
    startExpr_(rhs.startExpr_)
{}


Foam::expressions::exprResultStored::exprResultStored
(
    const dictionary& dict
)
:
    expressions::exprResult(dict.subOrEmptyDict("value")),
    name_(dict.get<word>("name")),
    startExpr_("initialValue", dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprResultStored::writeDict(Ostream& os) const
{
    os.beginBlock();

    os.writeEntry("name", name_);
    os.writeEntry("initialValue", startExpr_);

    os.writeKeyword("value");
    os << static_cast<const exprResult&>(*this);

    os.endBlock();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::expressions::exprResultStored::operator=
(
    const exprResultStored& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->exprResult::operator=(rhs);

    name_ = rhs.name_;
    startExpr_ = rhs.startExpr_;
}


void Foam::expressions::exprResultStored::operator=
(
    const exprResult& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->exprResult::operator=(rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    expressions::exprResultStored& data
)
{
    dictionary dict(is);
    data = expressions::exprResultStored(dict);

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprResultStored& data
)
{
    data.writeDict(os);

    return os;
}


// ************************************************************************* //
