/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "boundaryPatch.H"
#include "dictionary.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryPatch::boundaryPatch
(
    const word& name,
    const label index,
    const label size,
    const label start,
    const word& physicalType
)
:
    patchIdentifier(name, index, physicalType),
    size_(size),
    start_(start)
{}


Foam::boundaryPatch::boundaryPatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    patchIdentifier(name, dict, index),
    size_(dict.get<label>("nFaces")),
    start_(dict.get<label>("startFace"))
{}


Foam::boundaryPatch::boundaryPatch(const boundaryPatch& p, const label index)
:
    boundaryPatch(p)
{
    patchIdentifier::index() = index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryPatch::write(Ostream& os) const
{
    // if (!type_.empty()) os.writeEntry("type", type_);
    patchIdentifier::write(os);
    os.writeEntry("nFaces", size_);
    os.writeEntry("startFace", start_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundaryPatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
