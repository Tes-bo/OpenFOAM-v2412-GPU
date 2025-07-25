/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "IOmapDistribute.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(IOmapDistribute, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::IOmapDistribute::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOmapDistribute::IOmapDistribute(const IOobject& io)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistribute>();

    readIOcontents();
}


Foam::IOmapDistribute::IOmapDistribute
(
    const IOobject& io,
    const mapDistribute& map
)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistribute>();

    if (!readIOcontents())
    {
        mapDistribute::operator=(map);
    }
}


Foam::IOmapDistribute::IOmapDistribute
(
    const IOobject& io,
    mapDistribute&& map
)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistribute>();

    mapDistribute::transfer(map);

    readIOcontents();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOmapDistribute::readData(Istream& is)
{
    return (is >> *this).good();
}


bool Foam::IOmapDistribute::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// ************************************************************************* //
