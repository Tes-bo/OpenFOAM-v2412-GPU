/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "coordinateSystems.H"
#include "predicates.H"
#include "PtrListOps.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(coordinateSystems);
}

// File-local

//- Header name for 1806 and earlier
static const char* headerTypeCompat = "IOPtrList<coordinateSystem>";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::coordinateSystems::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // Attempt reading
    }
    else
    {
        return false;
    }


    // Do reading
    Istream& is = readStream(word::null);

    if (headerClassName() == typeName)
    {
        this->readIstream(is, coordinateSystem::iNew());
        close();
    }
    else if (headerClassName() == headerTypeCompat)
    {
        // Older (1806 and earlier) header name
        if (error::master())
        {
            std::cerr
                << "--> FOAM IOWarning :" << nl
                << "    Found header class name '" << headerTypeCompat
                << "' instead of '" << typeName << "'" << nl;

            error::warnAboutAge("header class", 1806);
        }

        this->readIstream(is, coordinateSystem::iNew());
        close();
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Unexpected class name " << headerClassName()
            << " expected " << typeName
            << " or " << headerTypeCompat << nl
            << "    while reading object " << name()
            << exit(FatalIOError);
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems(const IOobject& io)
:
    regIOobject(io)
{
    readIOcontents();
}


Foam::coordinateSystems::coordinateSystems(const objectRegistry& obr)
:
    coordinateSystems
    (
        IOobject
        (
            typeName,
            obr.time().constant(),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& content
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        static_cast<PtrList<coordinateSystem>&>(*this) = content;
    }
}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    PtrList<coordinateSystem>&& content
)
:
    regIOobject(io),
    PtrList<coordinateSystem>(std::move(content))
{
    readIOcontents();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::coordinateSystems& Foam::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    // Previously registered?

    const coordinateSystems* ptr = obr.findObject<coordinateSystems>(typeName);

    if (ptr)
    {
        return *ptr;
    }

    // Read construct from registry
    return obr.store(new coordinateSystems(obr));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::coordinateSystems::indices(const wordRe& key) const
{
    if (key.empty())
    {
        return labelList();
    }
    return PtrListOps::findMatching(*this, key);
}


Foam::labelList Foam::coordinateSystems::indices(const wordRes& matcher) const
{
    if (matcher.empty())
    {
        return labelList();
    }
    return PtrListOps::findMatching(*this, matcher);
}


Foam::label Foam::coordinateSystems::findIndex(const wordRe& key) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


Foam::label Foam::coordinateSystems::findIndex(const wordRes& matcher) const
{
    if (matcher.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, matcher);
}


bool Foam::coordinateSystems::found(const wordRe& key) const
{
    return findIndex(key) != -1;
}


const Foam::coordinateSystem*
Foam::coordinateSystems::cfind(const word& name) const
{
    const label index =
    (
        name.empty() ? -1 : PtrListOps::firstMatching(*this, name)
    );

    if (coordinateSystem::debug)
    {
        InfoInFunction
            << "Global coordinate system: "
            << name << '=' << index << endl;
    }

    // Return nullptr if not found
    return PtrList<coordinateSystem>::test(index);
}


const Foam::coordinateSystem&
Foam::coordinateSystems::lookup(const word& name) const
{
    const coordinateSystem* ptr = this->cfind(name);

    if (!ptr)
    {
        FatalErrorInFunction
            << "Could not find coordinate system: " << name << nl
            << "available coordinate systems: "
            << flatOutput(names()) << nl << nl
            << exit(FatalError);
    }

    return *ptr;
}


Foam::wordList Foam::coordinateSystems::names() const
{
    return PtrListOps::names(*this);  // match any/all
}


Foam::wordList Foam::coordinateSystems::names(const wordRe& key) const
{
    if (key.empty())
    {
        return wordList();
    }
    return PtrListOps::names(*this, key);
}


Foam::wordList Foam::coordinateSystems::names(const wordRes& matcher) const
{
    if (matcher.empty())
    {
        return wordList();
    }
    return PtrListOps::names(*this, matcher);
}


bool Foam::coordinateSystems::writeData(Ostream& os) const
{
    const PtrList<coordinateSystem>& list = *this;

    os << nl << size() << nl << token::BEGIN_LIST;

    for (const coordinateSystem& csys : list)
    {
        os << nl;
        csys.writeEntry(csys.name(), os);
    }

    os << token::END_LIST << nl;

    return os.good();
}


bool Foam::coordinateSystems::writeObject
(
    IOstreamOption,
    const bool writeOnProc
) const
{
    // Force ASCII, uncompressed
    return regIOobject::writeObject
    (
        IOstreamOption(IOstreamOption::ASCII),
        writeOnProc
    );
}


// ************************************************************************* //
