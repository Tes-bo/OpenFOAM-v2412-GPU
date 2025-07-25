/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::dictionary

Description
    A list of keyword definitions, which are a keyword followed by a number
    of values (eg, words and numbers) or by a sub-dictionary.
    Since the dictionary format is used extensively throughout OpenFOAM for
    input/output files, there are many examples of its use.

    Dictionary keywords are a plain word or a pattern (regular expression).
    The general order for searching is as follows:
      - exact match
      - pattern match (in reverse order)
      - optional recursion into the enclosing (parent) dictionaries

    The dictionary class is the base class for IOdictionary and also serves
    as a bootstrap dictionary for the objectRegistry data dictionaries.

Note
    Within dictionaries, entries can be referenced by using the
    \c '$' syntax familiar from shell programming.
    Similarly, the \c '/' separator is used when referencing
    sub-dictionary entries:
    - <b> "./" </b> : the current dictionary
    - <b> "../" </b> : the parent dictionary
    - <b> "../../" </b> : the grandparent dictionary
    .

    An initial \c '/' anchor specifies that the path starts from the
    top-level entry. It is also possible to use the '${}' syntax for clarity.

    For example,
    \verbatim
    key1        val1;
    key2        $key1;    // Use key1 value from current scope
    key3        $./key1;  // Use key1 value from current scope

    subdict1
    {
        key1        val1b;
        key2        $../key1;  // Use key1 value from parent
        subdict2
        {
            key2    val2;
            key3    $../../key1;  // Use key1 value from grandparent
        }
    }

    key4        $/subdict1/subdict2/key3;  // Lookup with absolute scoping
    \endverbatim

    Prior to OpenFOAM-v1712, a dot-scoping (.) syntax was used, which is
    still supported (AUG-2023) but deprecated in favour of the
    less ambiguous slash-scoping (/) syntax.
    With dot-scoping, an initial \c '^' anchor, or an initial (:), was used
    to specify that the path starts from the top-level entry.


SourceFiles
    dictionary.C
    dictionaryIO.C
    dictionarySearch.C

SeeAlso
    - Foam::entry
    - Foam::dictionaryEntry
    - Foam::primitiveEntry

\*---------------------------------------------------------------------------*/

#ifndef Foam_dictionary_H
#define Foam_dictionary_H

#include <type_traits>
#include "entry.H"
#include "IDLList.H"
#include "DLList.H"
#include "fileName.H"
#include "ITstream.H"
#include "HashTable.H"
#include "wordList.H"
#include "className.H"
#include "refPtr.H"
#include "IOobjectOption.H"

// Some common data types
#include "label.H"
#include "scalar.H"
#include "regExpFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class OSstream;
class SHA1Digest;

Istream& operator>>(Istream& is, dictionary& dict);
Ostream& operator<<(Ostream& os, const dictionary& dict);

/*---------------------------------------------------------------------------*\
                         Class dictionary Declaration
\*---------------------------------------------------------------------------*/

class dictionary
:
    public IDLList<entry>
{
public:

    // Searching

        //- Generic const/non-const dictionary entry %searcher.
        //  A %searcher provides a uniform means of finding and returning
        //  an entry pointer as well as the dictionary \a context in which
        //  the entry was located.
        //
        //  Note that the constructors and set methods are protected such
        //  that only friends of the class can set things. This safeguards
        //  against inconsistencies in context/entry.
        template<bool Const>
        class Searcher
        {
        public:
            friend dictionary;

            //- The const/non-const type for the context and sub-dictionaries
            typedef typename std::conditional
                <Const, const dictionary, dictionary>::type dict_type;

            //- The const/non-const type for entries
            typedef typename std::conditional
                <Const, const entry, entry>::type value_type;

            //- A pointer to a const/non-const dictionary
            typedef dict_type* dict_pointer;

            //- A reference to a const/non-const dictionary
            typedef dict_type& dict_reference;

            //- A pointer to a const/non-const entry
            typedef value_type* pointer;

            //- A reference to a const/non-const entry
            typedef value_type& reference;


        protected:

            //- The dictionary context for the entry
            dict_pointer dict_;

            //- The entry or nullptr
            pointer eptr_;


            //- Implicit construct for the given dictionary context
            Searcher(dict_pointer dict) noexcept
            :
                dict_(dict),
                eptr_(nullptr)
            {}

            //- Assign the entry
            void set(pointer eptr) noexcept
            {
                eptr_ = eptr;
            }


        public:

            //- Default construct
            constexpr Searcher() noexcept
            :
                dict_(nullptr),
                eptr_(nullptr)
            {}


            //- True if entry was found
            bool good() const noexcept { return eptr_; }

            //- True if entry was found
            //  \deprecated(2019-01) - prefer good() method
            bool found() const noexcept { return eptr_; }

            //- The containing dictionary context
            dict_reference context() const { return *dict_; }

            //- A pointer to the entry (nullptr if not found)
            pointer ptr() const noexcept { return eptr_; }

            //- A reference to the entry (Error if not found)
            reference ref() const { return *eptr_; }

            //- True if found entry is a dictionary.
            bool isDict() const noexcept
            {
                return (eptr_ && eptr_->dictPtr());
            }

            //- True if found entry is a stream.
            bool isStream() const noexcept
            {
                return (eptr_ && eptr_->streamPtr());
            }

            //- Pointer to the found entry as a dictionary, nullptr otherwise
            dict_pointer dictPtr() const noexcept
            {
                return (eptr_ ? eptr_->dictPtr() : nullptr);
            }

            //- Pointer to the found entry as a stream, nullptr otherwise
            ITstream* streamPtr() const noexcept
            {
                return (eptr_ ? eptr_->streamPtr() : nullptr);
            }

            //- Return the found entry as a dictionary.
            //- Error if not found, or not a dictionary.
            dict_reference dict() const { return eptr_->dict(); }

            //- Return the found entry as a ITstream.
            //- Error if not found, or not a stream.
            ITstream& stream() const { return eptr_->stream(); }

            //- Permit an explicit cast to the other (const/non-const) searcher
            explicit operator const Searcher<!Const>&() const
            {
                return *reinterpret_cast<const Searcher<!Const>*>(this);
            }

            //- A pointer to the entry (nullptr if not found)
            pointer operator->() const noexcept { return eptr_; }

            //- A reference to the entry (Error if not found)
            reference operator*() const { return *eptr_; }
        };


        //- Searcher with const access
        typedef Searcher<true> const_searcher;

        //- Searcher with non-const access
        typedef Searcher<false> searcher;


    // Friends

        //- Declare friendship with the entry class for IO
        friend class entry;

        //- Declare friendship with the searcher classes
        friend const_searcher;
        friend searcher;


private:

    // Private Data

        //- The dictionary name
        fileName name_;

        //- Parent dictionary
        const dictionary& parent_;

        //- Quick lookup of the entries held on the IDLList
        HashTable<entry*> hashedEntries_;

        //- Entries of matching patterns
        DLList<entry*> patterns_;

        //- Patterns as precompiled regular expressions
        DLList<autoPtr<regExp>> regexps_;


    // Typedefs

        //- The storage container
        typedef IDLList<entry> parent_type;


    // Private Member Functions

        //- Convert old-style (1806) boolean search specification to enum
        //
        //  \param recursive search parent dictionaries
        //  \param pattern match using regular expressions as well
        inline static enum keyType::option
        matchOpt(bool recursive, bool pattern)
        {
            return
                keyType::option
                (
                    (pattern ? keyType::REGEX : keyType::LITERAL)
                  | (recursive ? keyType::RECURSIVE : 0)
                );
        }

        //- Search using a '.' for scoping.
        //  A leading dot means to use the parent dictionary.
        //  An intermediate dot separates a sub-dictionary or sub-entry.
        //  However, the use of dots is unfortunately ambiguous.
        //  The value "a.b.c.d" could be a first-level entry, a second-level
        //  entry (eg, "a" with "b.c.d", "a.b" with "c.d" etc),
        //  or just about any other combination.
        //  The heuristic tries successively longer top-level entries
        //  until there is a suitable match.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt the search mode
        const_searcher csearchDotScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Search using a '/' for scoping.
        //  Semantics as per normal files: an intermediate "." is the current
        //  dictionary level, an intermediate ".." is the parent dictionary.
        //  Note that since a slash is not a valid word character, there is no
        //  ambiguity between separator and content.
        //  No possibility or need for recursion.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt the search mode. Recursive is ignored.
        const_searcher csearchSlashScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;


        //- Emit IOError about bad input for the entry
        void raiseBadInput(const ITstream& is, const word& keyword) const;

        //- The currently known executable name,
        //- obtained from argList envExecutable
        static word executableName();

        //- Report (usually stderr) that the keyword default value was used,
        //- or FatalIOError when writeOptionalEntries greater than 1
        template<class T>
        void reportDefault
        (
            const word& keyword,
            const T& deflt,
            const bool added = false  // Value was added to the dictionary
        ) const;


public:

    // Declare name of the class and its debug switch
    ClassName("dictionary");

    // Static Data

        //- Report optional keywords and values if not present in dictionary
        //  For value greater than 1: fatal.
        //  Set/unset via an InfoSwitch or -info-switch at the command-line
        static int writeOptionalEntries;

        //- An empty dictionary, which is also the parent for all dictionaries
        static const dictionary null;

        //- Output location when reporting default values
        static refPtr<OSstream> reportingOutput;


    // Static Member Functions

        //- Return the state of reporting optional (default) entries
        //  0: no reporting, 1: report, 2: fatal if not set
        inline static int reportOptional() noexcept;

        //- Change the state of reporting optional (default) entries
        //  0: no reporting, 1: report, 2: fatal if not set
        //  \return old level
        inline static int reportOptional(const int level) noexcept;


    // Constructors

        //- Default construct, a top-level empty dictionary
        dictionary();

        //- Construct top-level empty dictionary with given name
        explicit dictionary(const fileName& name);

        //- Construct given the entry name, parent dictionary and Istream,
        //- reading entries until EOF, optionally keeping the header
        dictionary
        (
            const fileName& name,
            const dictionary& parentDict,
            Istream& is,
            bool keepHeader = false
        );

        //- Construct top-level dictionary from Istream (discards the header).
        //- Reads entries until EOF or when the first token is a
        //- '{' character, it will stop reading at the matching '}' character.
        //  \note this constructor should be explicit
        dictionary(Istream& is);

        //- Construct top-level dictionary from Istream,
        //- reading entries until EOF, optionally keeping the header
        dictionary(Istream& is, bool keepHeader);

        //- Copy construct given the parent dictionary
        dictionary(const dictionary& parentDict, const dictionary& dict);

        //- Copy construct top-level dictionary
        dictionary(const dictionary& dict);

        //- Construct top-level dictionary as copy from pointer to dictionary.
        //  A null pointer is treated like an empty dictionary.
        explicit dictionary(const dictionary* dict);

        //- Move construct for given parent dictionary
        dictionary(const dictionary& parentDict, dictionary&& dict);

        //- Move construct top-level dictionary
        dictionary(dictionary&& dict);

        //- Construct and return clone
        autoPtr<dictionary> clone() const;

        //- Construct top-level dictionary on freestore from Istream
        static autoPtr<dictionary> New(Istream& is);


    //- Destructor
    virtual ~dictionary();


    // Member Functions

    // Access

        //- The dictionary name
        inline const fileName& name() const noexcept;

        //- The dictionary name for modification (use with caution).
        inline fileName& name() noexcept;

        //- The local dictionary name (final part of scoped name)
        inline word dictName() const;

        //- The dictionary name relative to the case.
        //  Uses argList::envRelativePath to obtain FOAM_CASE
        //
        //  \param caseTag replace globalPath with \<case\> for later
        //      use with expand(), or prefix \<case\> if the file name was
        //      not an absolute location
        fileName relativeName(const bool caseTag = false) const;

        //- The dictionary is actually dictionary::null (root dictionary)
        inline bool isNullDict() const noexcept;

        //- Return the parent dictionary
        inline const dictionary& parent() const noexcept;

        //- Return the top of the tree
        const dictionary& topDict() const;

        //- Return line number of first token in dictionary
        label startLineNumber() const;

        //- Return line number of last token in dictionary
        label endLineNumber() const;

        //- Return the SHA1 digest of the dictionary contents
        SHA1Digest digest() const;

        //- Return the dictionary as a list of tokens
        tokenList tokens() const;


    // Searching

        //- Find an entry (const access) with the given keyword.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to the entry found or a nullptr.
        inline const entry* findEntry
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find an entry (non-const access) with the given keyword.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to the entry found or a nullptr.
        inline entry* findEntry
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find an entry (const access) with the given keyword.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return True if entry was found
        inline bool found
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Search for a scoped entry (const access) with the given keyword.
        //  Allows scoping using '.'.
        //  Special handling for an absolute anchor (^) at start of the keyword
        //  and for '..' to ascend into the parent dictionaries.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to the entry found or a nullptr.
        inline const entry* findScoped
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a sub-dictionary pointer if present
        //- (and it is a dictionary) otherwise return nullptr.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to sub-dictionary found or a nullptr.
        inline const dictionary* findDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a sub-dictionary pointer if present
        //- (and it is a dictionary) otherwise return nullptr.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to sub-dictionary found or a nullptr.
        inline dictionary* findDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find and return an entry stream if present
        //- (and it is a stream) otherwise return nullptr.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return pointer to ITstream or a nullptr
        inline ITstream* findStream
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;


    // Lookup

        //- Search for an entry (const access) with the given keyword.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode
        //
        //  \return return an entry if present, otherwise FatalIOError.
        const entry& lookupEntry
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Find and return an entry data stream.
        //- FatalIOError if not found, or not a stream
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return entry stream or FatalIOError
        ITstream& lookup
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T.
        //- FatalIOError if not found, or if the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T>
        T get
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param deflt the default value to use
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T>
        T getOrDefault
        (
            const word& keyword,
            const T& deflt,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value
        //- and add it to dictionary.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param deflt the default value to use
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T>
        T getOrAdd
        (
            const word& keyword,
            const T& deflt,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find entry and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect,
        //- or it is mandatory and not found.
        //
        //  \param keyword the keyword to search for
        //  \param val the value to read into
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //  \param readOpt the entry is required/optional (default: MUST_READ)
        //
        //  \return true if the entry was read
        template<class T>
        bool readEntry
        (
            const word& keyword,
            T& val,
            enum keyType::option matchOpt = keyType::REGEX,
            IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ
        ) const;

        //- Find an entry if present, and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param val the value to read into
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return true if the entry was read
        template<class T>
        bool readIfPresent
        (
            const word& keyword,
            T& val,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T with additional checking
        //- FatalIOError if not found, or if the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param pred the value check predicate
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T, class Predicate>
        T getCheck
        (
            const word& keyword,
            const Predicate& pred,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param deflt the default value to use
        //  \param pred the value check predicate
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T, class Predicate>
        T getCheckOrDefault
        (
            const word& keyword,
            const T& deflt,
            const Predicate& pred,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value
        //- and add it to dictionary.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param deflt the default value to use
        //  \param pred the value check predicate
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T, class Predicate>
        T getCheckOrAdd
        (
            const word& keyword,
            const T& deflt,
            const Predicate& pred,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find entry and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect,
        //- or it is mandatory and not found.
        //
        //  \param keyword the keyword to search for
        //  \param val the value to read into
        //  \param pred the value check predicate
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //  \param readOpt the entry is required/optional (default: MUST_READ)
        //
        //  \return true if the entry was read
        template<class T, class Predicate>
        bool readCheck
        (
            const word& keyword,
            T& val,
            const Predicate& pred,
            enum keyType::option matchOpt = keyType::REGEX,
            IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ
        ) const;

        //- Find an entry if present, and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //
        //  \param keyword the keyword to search for
        //  \param val the value to read into
        //  \param pred the value check predicate
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return true if the entry read
        template<class T, class Predicate>
        bool readCheckIfPresent
        (
            const word& keyword,
            T& val,
            const Predicate& pred,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a sub-dictionary.
        //  Fatal if the entry does not exist or is not a sub-dictionary.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        const dictionary& subDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a sub-dictionary for manipulation.
        //  Fatal if the entry does not exist or is not a sub-dictionary.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        dictionary& subDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find and return a sub-dictionary for manipulation.
        //  Fatal if the entry exist and is not a sub-dictionary.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        dictionary& subDictOrAdd
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Find and return a sub-dictionary as a copy, otherwise return
        //- an empty dictionary.
        //  Warn if the entry exists but is not a sub-dictionary.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //  \param mandatory the keyword is mandatory (default: false)
        dictionary subOrEmptyDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX,
            const bool mandatory = false
        ) const;

        //- Find and return a sub-dictionary, otherwise return this dictionary.
        //  Warn if the entry exists but is not a sub-dictionary.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        const dictionary& optionalSubDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Return the table of contents
        wordList toc() const;

        //- Return the sorted table of contents
        wordList sortedToc() const;

        //- Return table of contents sorted using the specified comparator
        template<class Compare>
        wordList sortedToc(const Compare& comp) const;

        //- Return the list of available keys or patterns
        List<keyType> keys(bool patterns = false) const;


    // Editing

        //- Substitute the given keyword (which is prefixed by '$')
        //  with the corresponding sub-dictionary entries
        bool substituteKeyword
        (
            const word& keyword,
            bool mergeEntry = false
        );

        //- Substitute the given scoped keyword (which is prefixed by '$')
        //  with the corresponding sub-dictionary entries
        bool substituteScopedKeyword
        (
            const word& keyword,
            bool mergeEntry = false
        );

        //- Add a new entry.
        //
        //  \param entryPtr the entry to add
        //  \param mergeEntry dictionaries are interwoven and primitive
        //         entries are overwritten (default: false)
        //
        //  \return pointer to inserted entry, or place of merging
        //   or nullptr on failure
        entry* add(entry* entryPtr, bool mergeEntry=false);

        //- Add an entry.
        //
        //  \param e the entry to add
        //  \param mergeEntry dictionaries are interwoven and primitive
        //         entries are overwritten (default: false)
        //  \return pointer to inserted entry, or place of merging
        //   or nullptr on failure
        entry* add(const entry& e, bool mergeEntry=false);

        //- Add a word entry.
        //  \param overwrite force overwrite of an existing entry.
        //  \return pointer to inserted entry or nullptr on failure
        entry* add(const keyType& k, const word& v, bool overwrite=false);

        //- Add a string entry.
        //  \param overwrite force overwrite of an existing entry.
        //  \return pointer to inserted entry or nullptr on failure
        entry* add(const keyType& k, const string& v, bool overwrite=false);

        //- Add a label entry.
        //  \param overwrite force overwrite of an existing entry.
        //  \return pointer to inserted entry or nullptr on failure
        entry* add(const keyType& k, const label v, bool overwrite=false);

        //- Add a scalar entry.
        //  \param overwrite force overwrite of an existing entry.
        //  \return pointer to inserted entry or nullptr on failure
        entry* add(const keyType& k, const scalar v, bool overwrite=false);

        //- Add a dictionary entry.
        //  \param mergeEntry merge into an existing sub-dictionary
        //  \return pointer to inserted entry, or place of merging
        //   or nullptr on failure
        entry* add
        (
            const keyType& k,
            const dictionary& d,
            bool mergeEntry = false
        );

        //- Add a T entry
        //  \param overwrite force overwrite of existing entry
        //  \return pointer to inserted entry or nullptr on failure
        template<class T>
        entry* add(const keyType& k, const T& v, bool overwrite=false);

        //- Assign a new entry, overwriting any existing entry.
        //
        //  \return pointer to inserted entry or nullptr on failure
        entry* set(entry* entryPtr);

        //- Assign a new entry, overwriting any existing entry.
        //
        //  \return pointer to inserted entry or nullptr on failure
        entry* set(const entry& e);

        //- Assign a dictionary entry, overwriting any existing entry.
        //
        //  \return pointer to inserted entry or nullptr on failure
        entry* set(const keyType& k, const dictionary& v);

        //- Assign a T entry, overwriting any existing entry.
        //  \return pointer to inserted entry or nullptr on failure
        template<class T>
        entry* set(const keyType& k, const T& v);

        //- Remove an entry specified by keyword
        bool remove(const word& keyword);

        //- Change the keyword for an entry,
        //  \param overwrite force overwrite of an existing entry.
        bool changeKeyword
        (
            const keyType& oldKeyword,
            const keyType& newKeyword,
            bool overwrite=false
        );

        //- Merge entries from the given dictionary.
        //  Also merge sub-dictionaries as required.
        bool merge(const dictionary& dict);

        //- Clear the dictionary
        void clear();

        //- Transfer the contents of the argument and annul the argument.
        void transfer(dictionary& dict);


    // Read

        //- Check after reading if the input token stream has unconsumed
        //- tokens remaining or if there were no tokens in the first place.
        //  Emits FatalIOError
        void checkITstream(const ITstream& is, const word& keyword) const;

        //- Read dictionary from Istream (discards the header).
        //- Reads entries until EOF or when the first token is a
        //- '{' character, it will stop reading at the matching '}' character.
        bool read(Istream& is);

        //- Read dictionary from Istream (optionally keeping the header)
        //- Reads entries until EOF or when the first token is a
        //- '{' character, it will stop reading at the matching '}' character.
        bool read(Istream& is, bool keepHeader);


    // Write

        //- Write sub-dictionary with its dictName as its header
        void writeEntry(Ostream& os) const;

        //- Write sub-dictionary with the keyword as its header
        void writeEntry(const keyType& keyword, Ostream& os) const;

        //- Write dictionary entries.
        //  \param extraNewLine adds additional newline between entries
        //         for "top-level" dictionaries
        void writeEntries(Ostream& os, const bool extraNewLine=false) const;

        //- Write dictionary, normally with sub-dictionary formatting
        void write(Ostream& os, const bool subDict=true) const;


    // Searching

        //- Search dictionary for given keyword
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        const_searcher csearch
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Search dictionary for given keyword
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        const_searcher search
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Search dictionary for given keyword
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode (default: non-recursive with patterns)
        searcher search
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        );

        //- Search using scoping.
        //  There are two types of scoping available:
        //  -# dot-scoping, where a '.' is used to delineate the scope
        //  -# slash-scoping, where a '/' is used to delineate the scope
        //
        //  For dot-scoping, a leading '^' traverses to the top-level
        //  dictionary, leading dots mean use the parent dictionary and an
        //  intermediate dot separates a sub-dictionary or sub-entry.
        //  However, since the use of dots is ambiguous ("a.b.c" could be
        //  an entry itself or represent a "bc" entry from dictionary "a" etc),
        //  the heuristic backtracks and attempts successively longer
        //  top-level entries until a suitable match is found.
        //
        //  For slash-scoping, semantics similar to directory structures are
        //  used. A leading '/' traverses to the top-level dictionary,
        //  a single leading or intermediate '.' references the current
        //  dictionary level. A '..' pair references the parent dictionary.
        //  Any doubled slashes are silently ignored.
        //  Since a slash is not a valid keyword character, there is no
        //  ambiguity between separator and content.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode
        const_searcher csearchScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Search using dot or slash scoping.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode
        const_searcher searchScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        ) const;

        //- Search using dot or slash scoping.
        //
        //  \param keyword the keyword to search for
        //  \param matchOpt search mode
        searcher searchScoped
        (
            const word& keyword,
            enum keyType::option matchOpt
        );

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        const dictionary* cfindScopedDict(const fileName& dictPath) const;

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        const dictionary* findScopedDict(const fileName& dictPath) const;

        //- Locate a sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path does not exist
        dictionary* findScopedDict(const fileName& dictPath);

        //- Locate existing or create sub-dictionary using slash-scoping
        //  \return nullptr if the dictionary path could not be created
        dictionary* makeScopedDict(const fileName& dictPath);


    // Compatibility helpers

        //- Search dictionary for given keyword and any compatibility names
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //      Old version 1600=OpenFOAM-v3.0, 240=OpenFOAM-2.4.x,
        //      170=OpenFOAM-1.7.x,...
        //  \param matchOpt search mode (default: non-recursive with patterns)
        const_searcher csearchCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return an entry pointer if present, or return a nullptr,
        //- using any compatibility names if needed.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param matchOpt search mode
        const entry* findCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt
        ) const;

        //- Search dictionary for given keyword and any compatibility names
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param matchOpt search mode (default: non-recursive with patterns)
        bool foundCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return an entry if present, otherwise FatalIOError,
        //- using any compatibility names if needed.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param matchOpt search mode
        const entry& lookupEntryCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt
        ) const;

        //- Find and return an entry data stream,
        //- using any compatibility names if needed.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param matchOpt search mode (default: non-recursive with patterns)
        ITstream& lookupCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T
        //- using any compatibility names if needed.
        //- FatalIOError if not found, or if there are excess tokens.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T>
        T getCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find and return a T, or return the given default value
        //- using any compatibility names if needed.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param deflt the default value to use
        //  \param matchOpt search mode (default: non-recursive with patterns)
        template<class T>
        T getOrDefaultCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            const T& deflt,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;

        //- Find entry and assign to T val
        //- using any compatibility names if needed.
        //- FatalIOError if there are excess tokens.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param val the value to read
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //  \param readOpt the entry is required/optional (default: MUST_READ)
        //
        //  \return true if the entry was read
        template<class T>
        bool readCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            T& val,
            enum keyType::option matchOpt = keyType::REGEX,
            IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ
        ) const;

        //- Find an entry if present, and assign to T val
        //- using any compatibility names if needed.
        //- FatalIOError if it is found and there are excess tokens.
        //
        //  \param keyword the keyword to search for
        //  \param compat list of old compatibility keywords and the last
        //      OpenFOAM version for which they were used.
        //  \param val the value to read
        //  \param matchOpt search mode (default: non-recursive with patterns)
        //
        //  \return true if the entry was found.
        template<class T>
        bool readIfPresentCompat
        (
            const word& keyword,
            std::initializer_list<std::pair<const char*,int>> compat,
            T& val,
            enum keyType::option matchOpt = keyType::REGEX
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const dictionary& rhs);

        //- Include entries from the given dictionary.
        //  Warn, but do not overwrite existing entries.
        void operator+=(const dictionary& rhs);

        //- Conditionally include entries from the given dictionary.
        //  Do not overwrite existing entries.
        void operator|=(const dictionary& rhs);

        //- Unconditionally include entries from the given dictionary.
        //  Overwrite existing entries.
        void operator<<=(const dictionary& rhs);


    // IOstream operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream& is, dictionary& dict);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream& os, const dictionary& dict);


    // Shortcuts - when a templated classes also inherits from a dictionary

        #undef  defineDictionaryGetter
        #define defineDictionaryGetter(Func, Type)                            \
            /*! \brief Same as get\<Type\>(const word&, keyType::option) */   \
            Type Func                                                         \
            (                                                                 \
                const word& keyword,                                          \
                enum keyType::option matchOpt = keyType::REGEX                \
            ) const                                                           \
            {                                                                 \
                return get<Type>(keyword, matchOpt);                          \
            }

        defineDictionaryGetter(getBool, bool);
        defineDictionaryGetter(getLabel, label);
        defineDictionaryGetter(getScalar, scalar);
        defineDictionaryGetter(getString, string);
        defineDictionaryGetter(getWord, word);
        defineDictionaryGetter(getFileName, fileName);

        #undef defineDictionaryGetter


    // Housekeeping

        //- Check for existence of a sub-dictionary.
        //- Generally prefer findDict() for more flexibility.
        FOAM_DEPRECATED_STRICT(2024-05, "findDict()")
        bool isDict
        (
            const word& keyword,
            enum keyType::option matchOpt = keyType::REGEX
        ) const
        {
            return static_cast<bool>(findDict(keyword, matchOpt));
        }

        //- Same as getOrDefault()
        template<class T>
        FOAM_DEPRECATED_STRICT(2019-06, "getOrDefault()")
        T lookupOrDefault
        (
            const word& keyword,
            const T& deflt,
            enum keyType::option matchOpt = keyType::REGEX
        ) const
        {
            return getOrDefault<T>(keyword, deflt, matchOpt);
        }

        //- Same as getOrAdd()
        template<class T>
        FOAM_DEPRECATED_STRICT(2019-06, "getOrAdd()")
        T lookupOrAddDefault
        (
            const word& keyword,
            const T& deflt,
            enum keyType::option matchOpt = keyType::REGEX
        )
        {
            return getOrAdd<T>(keyword, deflt, matchOpt);
        }

        //- Deprecated(2018-07) - use lookup() method
        //  \deprecated(2018-07) - use lookup() method
        FOAM_DEPRECATED_FOR(2018-07, "lookup() method")
        ITstream& operator[](const word& keyword) const
        {
            return lookup(keyword);
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        FOAM_DEPRECATED_FOR(2018-10, "found(keyType::option)")
        bool found
        (
            const word& keyword,
            bool recursive,
            bool patternMatch = true
        ) const
        {
            return found(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use findEntry() method
        FOAM_DEPRECATED_FOR(2018-10, "findEntry(keyType::option)")
        entry* lookupEntryPtr
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        )
        {
            return findEntry(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use findEntry() method
        FOAM_DEPRECATED_FOR(2018-10, "findEntry(keyType::option)")
        const entry* lookupEntryPtr
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        ) const
        {
            return findEntry(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use findScoped() method
        FOAM_DEPRECATED_FOR(2018-10, "findScoped(keyType::option)")
        const entry* lookupScopedEntryPtr
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        ) const
        {
            return findScoped(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  Find and return a sub-dictionary pointer if present
        //  (and a sub-dictionary) otherwise return nullptr.
        //
        //  Search type: non-recursive with patterns.
        //  \deprecated(2018-10) - use findDict() method
        FOAM_DEPRECATED_FOR(2018-10, "findDict() method")
        const dictionary* subDictPtr(const word& keyword) const
        {
            return findDict(keyword, keyType::REGEX);
        }

        //- Deprecated(2018-10)
        //- Find and return a sub-dictionary pointer if present
        //  (and a sub-dictionary) otherwise return nullptr.
        //
        //  Search type: non-recursive with patterns.
        //  \deprecated(2018-10) - use findDict() method
        FOAM_DEPRECATED_FOR(2018-10, "findDict() method")
        dictionary* subDictPtr(const word& keyword)
        {
            return findDict(keyword, keyType::REGEX);
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        FOAM_DEPRECATED_FOR(2018-10, "lookupEntry(keyType::option)")
        const entry& lookupEntry
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        ) const
        {
            return lookupEntry(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        FOAM_DEPRECATED_FOR(2018-10, "lookup(keyType::option)")
        ITstream& lookup
        (
            const word& keyword,
            bool recursive,
            bool patternMatch = true
        ) const
        {
            return lookup(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        template<class T>
        FOAM_DEPRECATED_FOR(2018-10, "getOrDefault(keyType::option)")
        T lookupOrDefault
        (
            const word& keyword,
            const T& deflt,
            bool recursive,
            bool patternMatch = true
        ) const
        {
            return getOrDefault(keyword, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        template<class T>
        FOAM_DEPRECATED_FOR(2018-10, "getOrAdd(keyType::option)")
        T lookupOrAddDefault
        (
            const word& keyword,
            const T& deflt,
            bool recursive,
            bool patternMatch = true
        )
        {
            return getOrAdd(keyword, deflt, matchOpt(recursive, patternMatch));
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use keyType::option version
        template<class T>
        FOAM_DEPRECATED_FOR(2018-10, "readIfPresent(keyType::option)")
        bool readIfPresent
        (
            const word& keyword,
            T& val,
            bool recursive,
            bool patternMatch = true
        ) const
        {
            return
                readIfPresent
                (keyword, val, matchOpt(recursive, patternMatch));
        }


    // More compatibility

        //- Deprecated(2018-10) find and return a T.
        //  \deprecated(2018-10) - use get() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-10, "get() method")
        T lookupType
        (
            const word& keyword,
            bool recursive = false,
            bool patternMatch = true
        ) const
        {
            return get<T>(keyword, matchOpt(recursive, patternMatch));
        }

        #ifdef COMPAT_OPENFOAM_ORG
        //! Handle new openfoam-org method type (2019-11)
        //  Only available if compiled with COMPAT_OPENFOAM_ORG
        template<class T>
        FOAM_DEPRECATED_FOR(2019-11, "get() method - openfoam.org compat")
        T lookup
        (
            const word& keyword,
            bool recursive = false,
            bool patternMatch = true
        ) const
        {
            return get<T>(keyword, matchOpt(recursive, patternMatch));
        }
        #endif
};


// Global Operators

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Warn, but do not overwrite the entries from dict1.
dictionary operator+(const dictionary& dict1, const dictionary& dict2);

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Do not overwrite the entries from dict1.
dictionary operator|(const dictionary& dict1, const dictionary& dict2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionaryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dictionaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
