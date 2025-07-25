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
    Foam::entry

Description
    A keyword and a list of tokens is an 'entry'.

    An entry can be read, written and printed, and the types and values of
    its tokens analysed.  An entry is a high-level building block for data
    description.  It is a front-end for the token parser. A list of entries
    can be used as a set of keyword syntax elements, for example.

SourceFiles
    entry.C
    entryIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_entry_H
#define Foam_entry_H

#include "keyType.H"
#include "IDLList.H"
#include "fileName.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ITstream;
class dictionary;
class entry;
Ostream& operator<<(Ostream& os, const entry& e);


/*---------------------------------------------------------------------------*\
                           Class entry Declaration
\*---------------------------------------------------------------------------*/

class entry
:
    public IDLList<entry>::link
{
public:

    // Enums/Typedefs

        //- The input mode options
        enum class inputMode
        {
            MERGE,      //!< Merge sub-dictionaries when possible
            OVERWRITE,  //!< Keep last entry. Silently remove previous ones
            PROTECT,    //!< Keep initial entry. Silently ignore subsequent ones
            WARN,       //!< Keep initial entry. Warn about subsequent ones.
            ERROR,      //!< FatalError for duplicate entries
            GLOBAL      //!< Use global value from #globalInputMode variable
        };


private:

    // Private Data

        //- Keyword of entry
        keyType keyword_;


    // Private Member Functions

        //- Get the next valid keyword.
        //  \return True if it is a valid keyType.
        static bool getKeyword
        (
            keyType& keyword,
            token& keyToken,
            Istream& is
        );

        //- Get the next valid keyword.
        //  Warn when an invalid token is encountered that is not EOF or
        //  end-of-block
        //  \return True if it is a valid keyType.
        static bool getKeyword(keyType& keyword, Istream& is);

        //- Emit IOError about bad input for the entry
        void raiseBadInput(const ITstream& is) const;


protected:

    // Protected Member Functions

        //- Report a read warning (on std::cerr)
        static void reportReadWarning(const IOstream&, const std::string&);


public:

    //- Enable or disable use of function entries and variable expansions.
    static int disableFunctionEntries;

    //- The current global input-mode.
    static inputMode globalInputMode;


    // Constructors

        //- Construct from keyword
        entry(const keyType& keyword);

        //- Construct as copy
        entry(const entry& e);

        //- Construct on freestore as copy with reference to the
        //  dictionary the copy belongs to
        virtual autoPtr<entry> clone
        (
            const dictionary& parentDict
        ) const = 0;

        //- Construct on freestore as copy
        //  Note: the parent directory is set to dictionary::null
        virtual autoPtr<entry> clone() const;

        //- Construct from an Istream and insert into the dictionary
        //  \param parentDict dictionary to insert into
        //  \param is the input stream
        //  \param inpMode the input mode.
        //     The default is to use the currently active #globalInputMode
        //  \param endChar the expected end character (eg, a closing brace).
        //     The endChar is 0 if no expectations are asserted.
        static bool New
        (
            dictionary& parentDict,
            Istream& is,
            const inputMode inpMode = inputMode::GLOBAL,
            const int endChar = 0
        );

        //- Construct an entry from Istream.
        //  The expected input comprises a keyword followed by a
        //  dictionaryEntry or a primitiveEntry.
        //
        //  - The dictionaryEntry starts with a '{' left brace and ends
        //    with a '}' right brace.
        //  - The primitiveEntry ends with a ';' semi-colon.
        //
        // Example input
        // \verbatim
        //     key1 { ... }   // dictionary input
        //     key2 ... ;     // primitive input
        // \endverbatim
        //
        //  \return The #entry read, or nullptr on error.
        static autoPtr<entry> New(Istream& is);

        //- Reset the #globalInputMode to %merge
        static void resetInputMode();


    //- Destructor
    virtual ~entry() = default;


    // Member Functions

        //- Return keyword
        const keyType& keyword() const noexcept { return keyword_; }

        //- Return non-const access to keyword
        keyType& keyword() noexcept { return keyword_; }

        //- Return the entry name
        virtual const fileName& name() const = 0;

        //- Return the entry name for modification
        virtual fileName& name() = 0;

        //- Return the entry name relative to the current case
        virtual fileName relativeName() const = 0;

        //- Return line number of first token in dictionary
        virtual label startLineNumber() const = 0;

        //- Return line number of last token in dictionary
        virtual label endLineNumber() const = 0;


        //- True if this entry is a stream
        virtual bool isStream() const noexcept { return this->streamPtr(); }

        //- Return pointer to token stream, if it is a primitive entry,
        //- otherwise return nullptr
        virtual ITstream* streamPtr() const noexcept { return nullptr; }

        //- Return token stream, if entry is a primitive entry
        virtual ITstream& stream() const = 0;


        //- True if this entry is a dictionary
        virtual bool isDict() const noexcept { return this->dictPtr(); }

        //- Return pointer to dictionary, if entry is a dictionary,
        //- otherwise return nullptr.
        virtual const dictionary* dictPtr() const noexcept { return nullptr; }

        //- Return non-const pointer to dictionary, if entry is a dictionary,
        //- otherwise return nullptr.
        virtual dictionary* dictPtr() noexcept { return nullptr; }

        //- Return dictionary, if entry is a dictionary,
        //- otherwise Fatal.
        virtual const dictionary& dict() const = 0;

        //- Return non-const access to dictionary, if entry is a dictionary,
        //- otherwise Fatal.
        virtual dictionary& dict() = 0;


    // Read

        //- Check after reading if the input token stream has unconsumed
        //- tokens remaining or if there were no tokens in the first place.
        //  Emits FatalIOError
        void checkITstream(const ITstream& is) const;

        //- Get a T from the stream,
        //- FatalIOError if the number of tokens is incorrect.
        template<class T>
        T get() const
        {
            T val;
            readEntry<T>(val);
            return val;
        }

        //- Assign to T val,
        //- FatalIOError if the number of tokens is incorrect.
        //
        //  \param val the value to read into
        template<class T>
        void readEntry(T& val) const
        {
            ITstream& is = this->stream();
            is >> val;

            checkITstream(is);
        }

        //- Get a T from the stream,
        //- FatalIOError if the number of tokens is incorrect.
        //
        //  \param pred the value check predicate
        template<class T, class Predicate>
        T getCheck(const Predicate& pred) const
        {
            T val;
            readCheck<T>(val, pred);
            return val;
        }

        //- Assign to T val,
        //- FatalIOError if the number of tokens is incorrect.
        //
        //  \param val the value to read into
        //  \param pred the value check predicate
        template<class T, class Predicate>
        void readCheck(T& val, const Predicate& pred) const
        {
            ITstream& is = this->stream();
            is >> val;

            checkITstream(is);
            if (!pred(val))
            {
                raiseBadInput(is);
            }
        }


    // Write

        //- Write
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        void operator=(const entry& e);

        bool operator==(const entry& e) const;
        bool operator!=(const entry& e) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const entry& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
