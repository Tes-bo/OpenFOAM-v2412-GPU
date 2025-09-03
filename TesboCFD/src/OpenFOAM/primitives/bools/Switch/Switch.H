/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    true/false, on/off, yes/no, any/none.
    Also accepts 0/1 as a string and shortcuts t/f, y/n.

SourceFiles
    Switch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Switch_H
#define Foam_Switch_H

#include "bool.H"
#include "stdFoam.H"

// Avoid any pre-processor conflicts with enum names
#undef FALSE
#undef TRUE
#undef NO
#undef YES
#undef OFF
#undef ON
#undef NONE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class token;
class word;
class Switch;

// IOstream Operators

//- Read Switch from stream using Foam::Switch(Istream&)
Istream& operator>>(Istream& is, Switch& sw);

//- Write Switch to stream as its text value (eg, "true", "false")
Ostream& operator<<(Ostream& is, const Switch sw);

/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
public:

    // Data Types

        //- Switch enumerations corresponding to common text representations.
        //  \note The values here are critical for its proper behaviour.
        //  The values correspond to an index into the predefined output names
        //  for the c_str() method and the lower bit is tested for
        //  determining the true/false bool value.
        enum switchType : unsigned char
        {
            FALSE   = 0 /*!< "false" */, TRUE = 1 /*!< "true" */,
            NO      = 2 /*!< "no" */,    YES  = 3 /*!< "yes" */,
            OFF     = 4 /*!< "off" */,   ON   = 5 /*!< "on" */,
            NONE    = 6 /*!< "none" */,  ANY  = 7 /*!< "any" */,
            INVALID = 8 /*!< "invalid" (output only) */,
        };


private:

    // Private Data

        //- The logic and enumerated text representation stored in a byte
        unsigned char value_;


    // Private Member Functions

        //- Find switchType for given string. Return 'INVALID' if not found.
        //  With failOnError, trigger FatalError if not found
        static switchType parse(const std::string& str, const bool failOnError);


public:

    // Generated Methods

        //- Copy construct
        Switch(const Switch&) noexcept = default;

        //- Copy assignment
        Switch& operator=(const Switch&) noexcept = default;


    // Constructors

        //- Default construct as false
        constexpr Switch() noexcept
        :
            value_(switchType::FALSE)
        {}

        //- Implicit construct from enumerated value
        constexpr Switch(const switchType sw) noexcept
        :
            value_(sw)
        {}

        //- Implicit construct from bool
        constexpr Switch(const bool b) noexcept
        :
            value_(b ? switchType::TRUE : switchType::FALSE)
        {}

        //- Implicit construct from int (treat integer as bool value)
        constexpr Switch(const int i) noexcept
        :
            value_(i ? switchType::TRUE : switchType::FALSE)
        {}

        //- Construct from string - catches bad input.
        //  Use static find() method for a failsafe alternative
        explicit Switch(const std::string& str);

        //- Construct from character array - catches bad input.
        //  Use static find() method for a failsafe alternative
        explicit Switch(const char* str);

        //- Construct from float with rounding to zero given by
        //- the tolerance (default: 0.5)
        explicit Switch(const float val, const float tol=0.5);

        //- Construct from double with rounding to zero given by
        //- the tolerance (default: 0.5)
        explicit Switch(const double val, const double tol=0.5);

        //- Construct from token. Handles bool/label/word types.
        explicit Switch(const token& tok);

        //- Construct from dictionary lookup.
        //  FatalError if anything is incorrect.
        Switch
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict  //!< dictionary
        );

        //- Find the key in the dictionary and use the corresponding
        //- switch value or the default if not found in dictionary.
        //
        //  FatalIOError if the switch name is incorrect.
        //  Specifying warnOnly downgrades the FatalIOError to an IOWarning.
        Switch
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const Switch deflt,     //!< fallback if not found
            const bool warnOnly = false  //!< Warn (not fail) on bad input
        );

        //- Construct from Istream by reading a token
        explicit Switch(Istream& is);


    // Helpers

        //- Construct Switch from dictionary, with default value
        static Switch getOrDefault
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const Switch deflt = switchType::FALSE  //!< fallback if not found
        );

        //- Construct from dictionary, supplying default value so that if the
        //- value is not found, it is added into the dictionary.
        static Switch getOrAddToDict
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            dictionary& dict,       //!< dictionary
            const Switch deflt = switchType::FALSE //!< default value to add
        );


    // Static Member Functions

        //- A string representation of bool as "false" / "true"
        static const char* name(const bool b) noexcept;

        //- True if there is a switch type corresponding to the given string.
        static bool contains(const std::string& str);

        //- Find switchType for the given string, returning as a Switch that
        //- can be tested for good() or bad().
        static Switch find(const std::string& str);


    // Member Functions

        //- True if the Switch represents a valid enumeration
        bool good() const noexcept;

        //- True if the Switch does not represent a valid enumeration
        bool bad() const noexcept { return !good(); }

        //- The underlying enumeration value
        switchType type() const noexcept;

        //- Flip the type, so OFF becomes ON, etc.
        //  Ignored if the Switch is INVALID
        void negate() noexcept;

        //- A C-string representation of the Switch value
        const char* c_str() const noexcept;

        //- A string representation of the Switch value
        std::string str() const;

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict  //!< dictionary
        );


    // Member Operators

        //- Conversion to bool
        operator bool() const noexcept
        {
            return (value_ & 0x1);
        }

        //- Assignment from enumerated value
        Switch& operator=(const switchType sw) noexcept
        {
            value_ = sw;
            return *this;
        }

        //- Assignment from bool
        Switch& operator=(const bool b) noexcept
        {
            value_ = (b ? Switch::TRUE : Switch::FALSE);
            return *this;
        }


    // Housekeeping

        //- Same as contains()
        static bool found(const std::string& str) { return contains(str); }

        //- Deprecated(2020-01) From string with/without bad input test
        //  \deprecated(2020-01) - confusing syntax, use static find() method
        FOAM_DEPRECATED_FOR(2019-02, "static find() method")
        Switch(const std::string& str, bool allowBad);

        //- Deprecated(2020-01) From string with/without bad input test
        //  \deprecated(2020-01) - confusing syntax, use static find() method
        FOAM_DEPRECATED_FOR(2019-02, "static find() method")
        Switch(const char* str, bool allowBad);

        //- Deprecated(2020-01) Use good() method, or static contains() method
        //  \deprecated(2020-01) Use good() method, or static contains() method
        FOAM_DEPRECATED_FOR(2019-02, "good() or static contains() method")
        bool valid() const noexcept { return good(); }

        //- Same as getOrAddToDict()
        FOAM_DEPRECATED_STRICT(2019-06, "getOrAddToDict()")
        static Switch lookupOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const Switch deflt = switchType::FALSE
        )
        {
            return getOrAddToDict(name, dict, deflt);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
