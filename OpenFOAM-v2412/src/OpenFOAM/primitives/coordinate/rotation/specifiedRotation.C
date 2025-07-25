/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "specifiedRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

defineTypeName(specified);
//FUTURE addToRunTimeSelectionTable(coordinateRotation, specified, dictionary);

} // End namespace coordinateRotations
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::specified::specified()
:
    coordinateRotation(),
    Rmatrix_(sphericalTensor::I)
{}


Foam::coordinateRotations::specified::specified(const tensor& rot)
:
    coordinateRotation(),
    Rmatrix_(rot)
{}


Foam::coordinateRotations::specified::specified(const dictionary& dict)
:
    specified()
{
    // Not yet implemented, pending definition
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::specified::clear()
{
    Rmatrix_ = sphericalTensor::I;
}


Foam::tensor Foam::coordinateRotations::specified::R() const
{
    return Rmatrix_;
}


void Foam::coordinateRotations::specified::write(Ostream& os) const
{
    os << "specified rotation";
}


void Foam::coordinateRotations::specified::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());

    os.endBlock();
}


// ************************************************************************* //
