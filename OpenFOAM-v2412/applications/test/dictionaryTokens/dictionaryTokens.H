/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::dictionaryTokens

Description
    Provides a stream of tokens from a dictionary.

    This can be used to return a stream of tokens from a dictionary
    without overhead or needing to reparse information.

    For example,

    \code
    OPstream os = ...;
    dictionaryTokens toks(dict);

    while (toks.good())
    {
       os.write(*toks);
       ++toks;
    }

    \endcode
    Or alternatively,

    \code
    dictionaryTokens toks(dict);

    while (toks.good())
    {
       os << *toks << nl;
       ++toks;
    }
    \endcode

SourceFiles
    dictionaryTokens.C

\*---------------------------------------------------------------------------*/

#ifndef dictionaryTokens_H
#define dictionaryTokens_H

#include "dictionary.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "dictionaryListEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionaryTokens;


/*---------------------------------------------------------------------------*\
                      Class dictionaryTokens Declaration
\*---------------------------------------------------------------------------*/

class dictionaryTokens
{
public:

    // Forward declarations
    class primitive_iterator;
    class dictionary_iterator;


private:

    // Private Member Data

        //- Reference to the dictionary being streamed.
        const dictionary& dict_;

        //- The current entry iterator
        IDLList<entry>::const_iterator entryIter_;

        //- The current entry iterator for primitiveEntry types
        mutable autoPtr<primitive_iterator> primIter_;

        //- The current entry iterator for dictionaryEntry and
        //- dictionaryListEntry types
        mutable autoPtr<dictionary_iterator> dictIter_;

        //- Set/unset primitive and dictionary when changing to next entry
        bool setIterator() const;


    // Private Member Functions

        //- No copy construct
        dictionaryTokens(const dictionaryTokens&) = delete;

        //- No copy assignment
        void operator=(const dictionaryTokens&) = delete;


public:

    // Static Member Functions

        //- The entry keyword as word or string token
        static token keywordToken(const entry& e);


    // Constructors

        //- Construct from reference to dictionary to be streamed
        dictionaryTokens(const dictionary& dict);


    // Member Functions

        //- True if the token stream is in a valid state
        bool good() const;

        //- The current token, or undefined if the stream is in an invalid
        //- invalid state.
        const token& operator*() const;

        //- Advance to the next token and return the updated stream stream.
        bool operator++();

        //- Reset to beginning
        void rewind();

};


/*---------------------------------------------------------------------------*\
            Class dictionaryTokens::primitive_iterator Declaration
\*---------------------------------------------------------------------------*/

//- An iterator for a primitiveEntry
//
//  The token stream output has the form
//
//  \verbatim
//      keyword  content tokens ';'
//  \endverbatim
//
class dictionaryTokens::primitive_iterator
{
    // Private Member Data

        //- Reference to the tokenList being streamed.
        const tokenList* tokensPtr_;

        //- The keyword as a token (string, word or undefined)
        const token key_;

        //- The closing ';' as a token
        const token end_;

        //- The current position within the tokenList
        label pos_;


    // Private Member Functions

        //- No copy construct
        primitive_iterator(const primitive_iterator&) = delete;

        //- No copy assignment
        void operator=(const primitive_iterator&) = delete;

public:

    // Constructors

        //- Construct from reference to primitiveEntry
        primitive_iterator(const primitiveEntry& e);


    // Member Functions

        //- True if the entry has keyword or tokens and has not indexed beyond
        //- the final trailing ';'
        bool good() const;

        //- The current token, or undefined if the stream is invalid.
        const token& operator*() const;

        //- Advance to the next token and return the updated stream stream.
        bool operator++();

};


/*---------------------------------------------------------------------------*\
            Class dictionaryTokens::dictionary_iterator Declaration
\*---------------------------------------------------------------------------*/

//- An iterator for a dictionaryEntry and dictionaryListEntry
//
//  The token stream output has the form
//
//  \verbatim
//      keyword '{' content '}'
//  \endverbatim
//
//  or for the dictionaryListEntry the form
//
//  \verbatim
//      size '(' content ')'
//  \endverbatim
//
class dictionaryTokens::dictionary_iterator
{
    // Private Member Data

        //- The possible output states
        enum states { KEY=0, OPEN, CONTENT, CLOSE, END };

        //- The keyword or the size (dictionaryListEntry) as a token
        const token key_;

        //- The opening brace '{' or bracket '('
        const token lbrace_;

        //- The closing brace ')' or bracket ')'
        const token rbrace_;

        //- The current output state
        int state_;

        //- A streamer for the dictionary content
        dictionaryTokens dictTokens_;


    // Private Member Functions

        //- No copy construct
        dictionary_iterator(const dictionary_iterator&) = delete;

        //- No copy assignment
        void operator=(const dictionary_iterator&) = delete;


public:

    // Constructors

        //- Construct from reference to dictionaryEntry
        dictionary_iterator(const dictionaryEntry& e);

        //- Construct from reference to dictionaryListEntry
        dictionary_iterator(const dictionaryListEntry& e);


    // Member Functions

        //- In a valid state
        bool good() const;

        //- The current token, or undefined if the stream is invalid.
        const token& operator*() const;

        //- Advance to the next token and return the updated stream stream.
        bool operator++();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
