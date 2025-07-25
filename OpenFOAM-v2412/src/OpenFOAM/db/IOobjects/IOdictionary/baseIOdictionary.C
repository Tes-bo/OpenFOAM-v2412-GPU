/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "baseIOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(baseIOdictionary, 0);

    bool baseIOdictionary::writeDictionaries
    (
        debug::infoSwitch("writeDictionaries", 0)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::baseIOdictionary::baseIOdictionary
(
    const IOobject& io,
    const dictionary* fallback
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
}


Foam::baseIOdictionary::baseIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    baseIOdictionary(io, &dict)
{}


Foam::baseIOdictionary::baseIOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

const Foam::word& Foam::baseIOdictionary::name() const
{
    return regIOobject::name();
}


bool Foam::baseIOdictionary::readData(Istream& is)
{
    is >> *this;

    if (writeDictionaries && UPstream::master() && !is.bad())
    {
        Sout<< nl
            << "--- baseIOdictionary " << name()
            << ' ' << objectPath() << ":" << nl;
        writeHeader(Sout);
        writeData(Sout);
        Sout<< "--- End of baseIOdictionary " << name() << nl << endl;
    }

    return !is.bad();
}


bool Foam::baseIOdictionary::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::baseIOdictionary::operator=(const baseIOdictionary& rhs)
{
    dictionary::operator=(rhs);
}


// ************************************************************************* //
