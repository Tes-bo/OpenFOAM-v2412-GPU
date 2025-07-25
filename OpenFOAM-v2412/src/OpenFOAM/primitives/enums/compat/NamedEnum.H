/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::NamedEnum

Description
    A NamedEnum is a wrapper around a list of names that represent
    particular enumeration values.

    \deprecated(2017-05) This class is retained for compatibility only and
    should be NOT used for any new code.
    The Foam::Enum class is robuster, more flexible, easier to use.

See Also
    Foam::Enum

SourceFiles
    NamedEnum.C

\*---------------------------------------------------------------------------*/

#ifndef FoamCompat_NamedEnum_H
#define FoamCompat_NamedEnum_H

#include "HashTable.H"
#include "wordList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class dictionary;
template<class EnumType, int nEnum> class NamedEnum;

template<class EnumType, int nEnum>
Ostream& operator<<(Ostream& os, const NamedEnum<EnumType, nEnum>& wrapped);


/*---------------------------------------------------------------------------*\
                          Class NamedEnum Declaration
\*---------------------------------------------------------------------------*/

template<class EnumType, int nEnum>
class FOAM_DEPRECATED_FOR(2017-05, "Foam::Enum") NamedEnum
{
    //- The nEnum must be positive (non-zero)
    static_assert(nEnum > 0, "nEnum must be positive (non-zero)");


    // Private Member Data

        //- The values for the enum
        HashTable<int> lookup_;


    // Private Member Functions

        //- No copy construct
        NamedEnum(const NamedEnum&) = delete;

        //- No copy assignment
        void operator=(const NamedEnum&) = delete;


public:

    //- The type of enumeration wrapped by NamedEnum
    typedef EnumType value_type;


    // Static data members

        //- The set of names corresponding to the enumeration EnumType
        static const char* names[nEnum];


    // Constructors

        //- Construct from names
        NamedEnum();


    // Member Functions

    // Access

        //- The number of lookup names for the enumeration
        inline label size() const;

        //- The list of enum names
        inline wordList toc() const;

        //- The sorted list of enum names
        inline wordList sortedToc() const;

        //- The list of enum names, in construction order
        wordList words() const;

        //- The list of enum values, in construction order
        List<int> values() const;


    // Query

        //- Test if there is an enumeration corresponding to the given name.
        inline bool found(const word& enumName) const;

        //- Test if there is an enumeration corresponding to the given name.
        inline bool hasEnum(const word& enumName) const;

        //- Test if there is a name corresponding to the given enumeration.
        bool hasName(const EnumType e) const;


    // Lookup

        //- Lookup the key in the dictionary and return the corresponding
        // enumeration element based on its name.
        // Fatal if anything is incorrect.
        EnumType lookup
        (
            const word& key,
            const dictionary& dict
        ) const;

        //- Find the key in the dictionary and return the corresponding
        //  enumeration element based on its name.
        //  Return the default value if the key was not found in the dictionary.
        //  Fatal if enumerated name was incorrect.
        EnumType lookupOrDefault
        (
            const word& key,
            const dictionary& dict,
            const EnumType deflt
        ) const;


    // IO

        //- Read a word from Istream and return the corresponding enumeration
        EnumType read(Istream& is) const;

        //- Write the name representation of the enumeration to an Ostream
        //  A noop if the enumeration wasn't found.
        void write(const EnumType e, Ostream& os) const;


    // Member Operators

        //- Return the enumeration element corresponding to the given name
        inline const EnumType operator[](const word& name) const;

        //- Return the name of the given enumeration element
        inline const char* operator[](const EnumType e) const;


    // IOstream operators

        //- Write names to Ostream, as per writeKeys() with shortListLen=10
        friend Ostream& operator<< <EnumType, nEnum>
        (
            Ostream& os,
            const NamedEnum<EnumType, nEnum>& wrapped
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NamedEnumI.H"

#ifdef NoRepository
    #include "NamedEnum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
