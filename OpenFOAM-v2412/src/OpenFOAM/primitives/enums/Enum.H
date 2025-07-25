/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::Enum

Description
    Enum is a wrapper around a list of names/values that represent particular
    enumeration (or int) values.
    All dictionary searches use a literal (not regex).

SourceFiles
    Enum.C
    EnumI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Enum_H
#define Foam_Enum_H

#include "wordList.H"
#include <ostream>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
template<class EnumType> class Enum;

/*---------------------------------------------------------------------------*\
                            Class Enum Declaration
\*---------------------------------------------------------------------------*/

template<class EnumType>
class Enum
{
    // Private Member Data

        //- The names for the enum
        List<word> keys_;

        //- The values for the enum, stored as int
        List<int> vals_;


    // Allow enums and integrals (should fit within an int)
    static_assert
    (
        std::is_enum<EnumType>::value || std::is_integral<EnumType>::value,
        "Enum must be enum or an integral type"
    );

public:

    // Typedefs

        //- The type of keys used
        typedef word key_type;

        //- The type of enumeration represented by the Enum
        typedef EnumType value_type;


    // Constructors

        //- Default construct, an empty list
        Enum() noexcept = default;

        //- Construct from a values/names list.
        //  Duplicate values are permitted (eg, for aliases).
        //  Duplicate names are permitted, but won't make much sense.
        explicit Enum
        (
            std::initializer_list<std::pair<EnumType, const char*>> list
        );


    // Member Functions

    // Access

        //- True if the enumeration list is empty.
        inline bool empty() const noexcept;

        //- The number of name/value pairs for the enumeration.
        inline label size() const noexcept;

        //- The list of enum names, in construction order. Same as toc()
        inline const List<word>& names() const noexcept;

        //- The list of enum values, in construction order.
        inline const List<int>& values() const noexcept;

        //- The list of enum names, in construction order. Same as names()
        inline const List<word>& toc() const noexcept;

        //- The sorted list of enum names.
        inline List<word> sortedToc() const;


    // Modify

        //- Clear all entries
        inline void clear();

        //- Append value/key pairs to the lists of known enumerations
        //  Does not check for duplicate entries
        void push_back
        (
            std::initializer_list<std::pair<EnumType, const char*>> list
        );


    // Query

        //- True if there is an enumeration corresponding to the given name.
        inline bool contains(const word& enumName) const;

        //- True if there is a name corresponding to the given enumeration.
        inline bool contains(const EnumType e) const;

        //- The enumeration corresponding to the given name.
        //  FatalError if not found.
        EnumType get(const word& enumName) const;

        //- The name corresponding to the given enumeration.
        //  Return an empty word if there is no corresponding name for it.
        inline const word& get(const EnumType e) const;

        //- The enumeration corresponding to the given name.
        //  \return The enumeration or default if not found.
        //  \note Method name compatibility with HashTable
        EnumType lookup(const word& enumName, const EnumType deflt) const;


    // Read

        //- Get the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //  FatalIOError if anything is incorrect.
        EnumType get
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict  //!< dictionary
        ) const;

        //- Find the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //
        //  \return The value found or default if not found in dictionary.
        //  FatalIOError if the enumeration is incorrect.
        //  Specifying warnOnly downgrades the FatalIOError to an IOWarning.
        EnumType getOrDefault
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const EnumType deflt,   //!< fallback if not found
            const bool warnOnly = false  //!< Warn (not fail) on bad enumeration
        ) const;

        //- Find entry and assign to T val.
        //  FatalIOError if the enumeration is incorrect,
        //  or when it is mandatory but was not found.
        //
        //  (TDB: handle IOobjectOption::readOption instead)
        //  \return true if the entry was found.
        bool readEntry
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            EnumType& val,          //!< The value to read into
            const bool mandatory = true, //!< The keyword is mandatory
            const bool warnOnly = false  //!< Warn (not fail) on bad enumeration
        ) const;

        //- Find an entry if present, and assign to T val.
        //  FatalIOError if the enumeration is incorrect.
        //  Default search: non-recursive with patterns.
        //
        //  \return true if the entry was found.
        inline bool readIfPresent
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            EnumType& val,          //!< The value to read into
            const bool warnOnly = false  //!< Warn (not fail) on bad enumeration
        ) const;


    // IO

        //- Read a word from Istream and return the corresponding enumeration
        EnumType read(Istream& is) const;

        //- Read a word from Istream, lookup named enumeration.
        //  \return true on success. Fatal if mandatory and not found.
        bool read
        (
            Istream& is,
            EnumType& val,
            const bool mandatory = true
        ) const;

        //- Write the name representation of the enumeration to an Ostream
        //  A no-op if the enumeration does not have a corresponding name.
        inline void write(const EnumType e, Ostream& os) const;

        //- Write enumeration names as a list without line-breaks
        //- to an output stream.
        template<class OS>
        inline OS& writeList(OS& os, const label ununsed=0) const;


    // Member Operators

        //- Return the enumeration corresponding to the given name
        //  FatalError if the name is not found.
        //  Identical to get()
        inline EnumType operator[](const word& enumName) const;

        //- Return the first name corresponding to the given enumeration,
        //- or an empty word on failure.
        //  Identical to get()
        inline const word& operator[](const EnumType e) const;


    // Iteration

        //- A const_iterator for iterating an Enum list
        //  \note The iterator dereference returns the \b key
        class const_iterator
        {
            //- The list being iterated
            const Enum* ptr_;

            //- Index in the list
            label idx_;

        public:

            //- Default construct, construct at given position
            inline explicit const_iterator
            (
                const Enum* eptr = nullptr,
                const label idx = 0
            ) noexcept;

            //- True if iterator points to an entry
            inline bool good() const noexcept;

            //- The name at the current index
            inline const word& key() const;

            //- Enumeration value at the current index
            inline EnumType val() const;

            //- De-referencing returns the name (key)
            //  This is similar to HashSet (not HashTable!) and allows
            //  convenient output and traversing of the names
            const word& operator*() const { return key(); }

            //- Move to the next index
            inline const_iterator& operator++() noexcept;

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };

        inline const_iterator cbegin() const noexcept;
        inline const_iterator cend() const noexcept;

        const_iterator begin() const noexcept { return cbegin(); }
        const_iterator end() const noexcept { return cend(); }

        //- Find key/value pair by enumeration name.
        //  Equal to cend() if not found, or test if good()
        inline const_iterator cfind(const word& key) const;
        const_iterator find(const word& key) const { return cfind(key); }

        //- Find key/value pair by enumeration value.
        //  Equal to cend() if not found, or test if good()
        inline const_iterator cfind(const EnumType e) const;
        const_iterator find(const EnumType e) const { return cfind(e); }


    // Housekeeping

        //- Same as getOrDefault()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrDefault() method")
        EnumType lookupOrDefault
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const EnumType deflt,   //!< fallback if not found
            const bool warnOnly = false  //!< Warn (not fail) on bad enumeration
        ) const
        {
            return getOrDefault(key, dict, deflt, warnOnly);
        }

        //- Deprecated(2020-11) use get() method
        //  \deprecated(2020-11) - use get() method
        FOAM_DEPRECATED_FOR(2020-11, "get() method")
        const word& operator()(const EnumType e) const
        {
            return get(e);
        }

        //- Deprecated(2020-11) use two-parameter lookup() method
        //  \deprecated(2020-11) - use two-parameter lookup() method
        FOAM_DEPRECATED_FOR(2020-11, "lookup() method")
        EnumType operator()(const word& key, const EnumType deflt) const
        {
            return lookup(key, deflt);
        }

        //- Deprecated(2020-11) use two-parameter lookup() method
        //  \deprecated(2020-11) - use two-parameter lookup() method
        FOAM_DEPRECATED_FOR(2020-11, "lookup() method")
        EnumType get(const word& key, const EnumType deflt) const
        {
            return lookup(key, deflt);
        }

        //- Deprecated(2018-10) same as two-parameter get() method
        //  \deprecated(2018-10) - use two-parameter get() method
        FOAM_DEPRECATED_FOR(2018-10, "get() method")
        EnumType lookup(const word& key, const dictionary& dict) const
        {
            return get(key, dict);
        }

        //- Append value/key pairs to the lists of known enumerations
        //  Does not check for duplicate entries
        FOAM_DEPRECATED_STRICT(2023-06, "push_back() method")
        void append
        (
            std::initializer_list<std::pair<EnumType, const char*>> list
        )
        {
            push_back(list);
        }

        //- Same as contains()
        bool found(const word& key) const { return contains(key); }

        //- Same as contains()
        bool found(const EnumType e) const { return contains(e); }
};


// Ostream Operator

//- Write enumeration names, without line-breaks (ie, FlatOutput)
template<class EnumType>
inline Ostream& operator<<(Ostream& os, const Enum<EnumType>& list);

//- Write enumeration names, without line-breaks (ie, FlatOutput)
template<class EnumType>
inline std::ostream& operator<<(std::ostream& os, const Enum<EnumType>& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EnumI.H"

#ifdef NoRepository
    #include "Enum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
