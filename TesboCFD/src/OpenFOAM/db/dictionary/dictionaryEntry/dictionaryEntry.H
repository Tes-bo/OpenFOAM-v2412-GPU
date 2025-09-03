/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

Class
    Foam::dictionaryEntry

Description
    A keyword and a list of tokens is a 'dictionaryEntry'.

    An dictionaryEntry can be read, written and printed, and the types and
    values of its tokens analysed.  A dictionaryEntry is a high-level building
    block for data description.  It is a front-end for the token parser.
    A list of entries can be used as a set of keyword syntax elements,
    for example.

SourceFiles
    dictionaryEntry.C
    dictionaryEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dictionaryEntry_H
#define Foam_dictionaryEntry_H

#include "entry.H"
#include "dictionary.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionaryEntry;
Ostream& operator<<(Ostream& os, const dictionaryEntry& e);

template<>
Ostream& operator<<(Ostream&, const InfoProxy<dictionaryEntry>&);


/*---------------------------------------------------------------------------*\
                       Class dictionaryEntry Declaration
\*---------------------------------------------------------------------------*/

class dictionaryEntry
:
    public entry,
    public dictionary
{
public:

    // Generated Methods

        //- No copy construct
        dictionaryEntry(const dictionaryEntry&) = delete;


    // Constructors

        //- Construct from the parent dictionary and Istream.
        //  The keyword is extracted from the stream
        dictionaryEntry(const dictionary& parentDict, Istream& is);

        //- Construct from the keyword, parent dictionary and a Istream
        dictionaryEntry
        (
            const keyType& key,
            const dictionary& parentDict,
            Istream& is
        );

        //- Construct from the keyword, parent dictionary and a dictionary
        dictionaryEntry
        (
            const keyType& key,
            const dictionary& parentDict,
            const dictionary& dict
        );

        //- Construct as copy for the given parent dictionary
        dictionaryEntry
        (
            const dictionary& parentDict,
            const dictionaryEntry& dictEnt
        );

        autoPtr<entry> clone(const dictionary& parentDict) const
        {
            return autoPtr<entry>(new dictionaryEntry(parentDict, *this));
        }


    // Member functions

        //- Return the scoped dictionary name (eg, dictA.dictB.dictC)
        virtual const fileName& name() const
        {
            return dictionary::name();
        }

        //- Return scoped dictionary name for modification
        virtual fileName& name()
        {
            return dictionary::name();
        }

        //- Return scoped dictionary name relative to the current case
        virtual fileName relativeName() const
        {
            return dictionary::relativeName();
        }

        //- Return line number of first token in dictionary
        virtual label startLineNumber() const;

        //- Return line number of last token in dictionary
        virtual label endLineNumber() const;

        //- This entry is not a primitive,
        //- calling this function generates a FatalError
        virtual ITstream& stream() const;


        //- Return pointer to this dictionary
        virtual const dictionary* dictPtr() const noexcept;

        //- Return non-const pointer to this dictionary
        virtual dictionary* dictPtr() noexcept;

        //- Return dictionary (ie, this)
        virtual const dictionary& dict() const noexcept;

        //- Return non-const access to dictionary
        virtual dictionary& dict() noexcept;


        //- Write
        virtual void write(Ostream& os) const;

        //- Return info proxy,
        //- used to print token information to a stream
        InfoProxy<dictionaryEntry> info() const noexcept { return *this; }


    // Ostream Operators

        friend Ostream& operator<<(Ostream& os, const dictionaryEntry& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
