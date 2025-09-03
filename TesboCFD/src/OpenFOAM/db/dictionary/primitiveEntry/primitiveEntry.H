/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::primitiveEntry

Description
    A keyword and a list of tokens comprise a primitiveEntry.
    A primitiveEntry can be read, written and printed, and the types and
    values of its tokens analysed.

    A primitiveEntry is a high-level building block for data description.
    It is a front-end for the token parser. A list of entries can be used
    as a set of keyword syntax elements, for example.

SourceFiles
    primitiveEntry.C
    primitiveEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitiveEntry_H
#define Foam_primitiveEntry_H

#include "entry.H"
#include "ITstream.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class primitiveEntry;

template<>
Ostream& operator<<(Ostream&, const InfoProxy<primitiveEntry>&);


/*---------------------------------------------------------------------------*\
                       Class primitiveEntry Declaration
\*---------------------------------------------------------------------------*/

class primitiveEntry
:
    public entry,
    public ITstream
{
    // Private Member Functions

        //- Test if token is acceptable after filtering for function entries
        //- and variable expansions.
        bool acceptToken
        (
            const token& tok,
            const dictionary& dict,
            Istream& is
        );

        //- Expand the given variable.
        //  The keyword starts with '$', but has been removed by the caller
        //  and thus passed as a varName.
        //  Keywords with '${}' are expanded recursively.
        bool expandVariable
        (
            const string& varName,
            const dictionary& dict
        );

        //- Expand the given function.
        //  The keyword starts with '#', but has been removed by the caller.
        //  and thus passed as a functionName.
        bool expandFunction
        (
            const word& functionName,
            const dictionary& dict,
            Istream& is
        );

        //- Read the complete entry from the given stream
        void readEntry(const dictionary& dict, Istream& is);


public:

    // Constructors

        //- Construct from keyword and no tokens.
        //  Contents to be filled with a later assignment
        explicit primitiveEntry(const keyType& key);

        //- Construct from keyword and a single token (copy)
        primitiveEntry(const keyType& key, const token& tok);

        //- Construct from keyword and a single token (move)
        primitiveEntry(const keyType& key, token&& tok);

        //- Construct from keyword and a list of tokens
        primitiveEntry(const keyType& key, const UList<token>& tokens);

        //- Construct from keyword and by transferring a list of tokens
        primitiveEntry(const keyType& key, List<token>&& tokens);

        //- Construct from keyword and ITstream tokens
        primitiveEntry(const keyType& key, const ITstream& is);

        //- Construct from keyword and Istream
        primitiveEntry(const keyType& key, Istream& is);

        //- Construct from keyword, parent dictionary and Istream
        primitiveEntry(const keyType& key, const dictionary& dict, Istream& is);

        //- Construct from keyword and a value. Uses string stream serialization
        template<class T>
        primitiveEntry(const keyType& key, const T& val);

        //- Clone the entry
        autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new primitiveEntry(*this));
        }


    // Member Functions

        //- Inherit read from ITstream
        using ITstream::read;

        //- Return the token stream name
        virtual const fileName& name() const
        {
            return ITstream::name();
        }

        //- Return token stream name for modification
        virtual fileName& name()
        {
            return ITstream::name();
        }

        //- Return token stream name relative to the current case
        virtual fileName relativeName() const
        {
            return ITstream::relativeName();
        }

        //- The line number of the first token in the entry
        virtual label startLineNumber() const
        {
            return ITstream::startLineNumber();
        }

        //- The line number of the last token in the entry
        virtual label endLineNumber() const
        {
            return ITstream::endLineNumber();
        }

        //- Return pointer to token stream for this primitive entry
        virtual ITstream* streamPtr() const noexcept;

        //- Return token stream for this primitive entry
        virtual ITstream& stream() const;

        //- This entry is not a dictionary,
        //- calling this function generates a FatalError
        virtual const dictionary& dict() const;

        //- This entry is not a dictionary,
        //- calling this function generates a FatalError
        virtual dictionary& dict();

        //- Read tokens from the given stream
        virtual bool read(const dictionary& dict, Istream& is);

        //- Write
        virtual void write(Ostream& os) const;

        //- Write, optionally with contents only (no keyword, etc)
        void write(Ostream& os, const bool contentsOnly) const;

        //- Return info proxy,
        //- to print token information to a stream
        InfoProxy<primitiveEntry> info() const noexcept { return *this; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitiveEntryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
