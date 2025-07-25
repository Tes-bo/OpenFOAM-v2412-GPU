/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    const dictionary* fallback
)
:
    IOdictionary(io, typeName, fallback)
{}


Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    IOdictionary(io, typeName, &dict)
{}


Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    const word& wantedType,
    const dictionary* fallback
)
:
    baseIOdictionary(io, fallback)
{
    if (!readHeaderOk(IOstreamOption::ASCII, wantedType) && fallback)
    {
        dictionary::operator=(*fallback);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    baseIOdictionary(io, is)
{
    // Default construct dictionary and read in afterwards
    // so that if there is some fancy massaging due to a
    // functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::dictionary Foam::IOdictionary::readContents(const IOobject& io)
{
    return readContents(io, typeName);
}


Foam::dictionary Foam::IOdictionary::readContents
(
    const IOobject& io,
    const word& wantedType
)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    // The object is global
    rio.globalObject(true);

    IOdictionary reader
    (
        rio,
        (wantedType.empty() ? typeName : wantedType)
    );

    return dictionary(std::move(static_cast<dictionary&>(reader)));
}


// ************************************************************************* //
