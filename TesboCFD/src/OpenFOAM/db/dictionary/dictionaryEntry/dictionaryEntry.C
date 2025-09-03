/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    const dictionary& dict
)
:
    entry(key),
    dictionary(parentDict, dict)
{}


Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    const dictionaryEntry& dictEnt
)
:
    entry(dictEnt),
    dictionary(parentDict, dictEnt)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dictionaryEntry::startLineNumber() const
{
    return dictionary::startLineNumber();
}


Foam::label Foam::dictionaryEntry::endLineNumber() const
{
    return dictionary::endLineNumber();
}


Foam::ITstream& Foam::dictionaryEntry::stream() const
{
    FatalIOErrorInFunction(*this)
        << "Attempt to return stream of primitives from a dictionary entry: "
        << entry::keyword() << nl
        << abort(FatalIOError);

    // Need to return something - send back an empty stream
    return ITstream::empty_stream();
}


const Foam::dictionary* Foam::dictionaryEntry::dictPtr() const noexcept
{
    return this;
}


Foam::dictionary* Foam::dictionaryEntry::dictPtr() noexcept
{
    return this;
}


const Foam::dictionary& Foam::dictionaryEntry::dict() const noexcept
{
    return *this;
}


Foam::dictionary& Foam::dictionaryEntry::dict() noexcept
{
    return *this;
}


// ************************************************************************* //
