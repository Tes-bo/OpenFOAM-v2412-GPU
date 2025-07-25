/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "STARCDCoordinateRotation.H"
#include "EulerCoordinateRotation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

    defineTypeName(starcd);

    // Standard short name
    addNamedToRunTimeSelectionTable
    (
        coordinateRotation,
        starcd,
        dictionary,
        starcd
    );

    // Long name - Compat 1806
    addAliasToRunTimeSelectionTable
    (
        coordinateRotation,
        starcd,
        dictionary,
        starcd,
        STARCDRotation,
        1806
    );

} // End namespace coordinateRotation
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::starcd::rotation
(
    const vector& angles,
    bool degrees
)
{
    return euler::rotation(euler::eulerOrder::ZXY, angles, degrees);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::starcd::starcd()
:
    coordinateRotation(),
    angles_(Zero),
    degrees_(true)
{}


Foam::coordinateRotations::starcd::starcd(const starcd& crot)
:
    coordinateRotation(),
    angles_(crot.angles_),
    degrees_(crot.degrees_)
{}


Foam::coordinateRotations::starcd::starcd
(
    const vector& rotZrotXrotY,
    bool degrees
)
:
    coordinateRotation(),
    angles_(rotZrotXrotY),
    degrees_(degrees)
{}


Foam::coordinateRotations::starcd::starcd
(
    scalar rotZ,
    scalar rotX,
    scalar rotY,
    bool degrees
)
:
    coordinateRotation(),
    angles_(rotZ, rotX, rotY),
    degrees_(degrees)
{}


Foam::coordinateRotations::starcd::starcd(const dictionary& dict)
:
    coordinateRotation(),
    angles_(dict.get<vector>("angles")),
    degrees_(dict.getOrDefault("degrees", true))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::starcd::clear()
{
    angles_ = Zero;
    degrees_ = true;
}


Foam::tensor Foam::coordinateRotations::starcd::R() const
{
    return starcd::rotation(angles_, degrees_);
}


void Foam::coordinateRotations::starcd::write(Ostream& os) const
{
    os  << "starcd-angles(" << (degrees_ ? "deg" : "rad") << "): " << angles_;
}


void Foam::coordinateRotations::starcd::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());
    os.writeEntry("angles", angles_);
    if (!degrees_)
    {
        os.writeEntry("degrees", "false");
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
