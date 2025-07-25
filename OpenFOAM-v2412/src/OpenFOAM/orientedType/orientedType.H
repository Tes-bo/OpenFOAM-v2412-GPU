/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::orientedType

Description
    Class to determine the 'oriented' status of surface fields

SourceFiles
    orientedType.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_orientedType_H
#define Foam_orientedType_H

#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class orientedType;

Istream& operator>>(Istream& is, orientedType& ot);
Ostream& operator<<(Ostream& os, const orientedType& ot);

/*---------------------------------------------------------------------------*\
                        Class orientedType Declaration
\*---------------------------------------------------------------------------*/

class orientedType
{
public:

    // Public Data Types

        //- Enumeration defining oriented flags
        enum orientedOption : unsigned char
        {
            UNKNOWN = 0,        //!< Unknown/undefined orientation
            ORIENTED = 1,       //!< Is oriented
            UNORIENTED = 2      //!< Is unoriented
        };


    // Static Data

        //- Named enumerations for oriented flags
        static const Enum<orientedOption> orientedOptionNames;


private:

    // Private Data

        //- Oriented type
        orientedOption oriented_;


public:

    // Generated Methods

        //- Copy construct
        orientedType(const orientedType&) noexcept = default;

        //- Copy assignment
        orientedType& operator=(const orientedType&) noexcept = default;


    // Constructors

        //- Default construct as \c UNKNOWN
        constexpr orientedType() noexcept
        :
            oriented_(UNKNOWN)
        {}

        //- Implicit construct from enumeration
        orientedType(orientedOption opt) noexcept
        :
            oriented_(opt)
        {}

        //- Construct from bool
        explicit constexpr orientedType(const bool isOriented) noexcept
        :
            oriented_(isOriented ? ORIENTED : UNORIENTED)
        {}

        //- Construct from Istream
        explicit orientedType(Istream& is);


    // Member Functions

        //- True if can operate on this combination of oriented types
        static bool checkType
        (
            const orientedType& a,
            const orientedType& b
        ) noexcept;

        //- Return the oriented flag
        orientedOption oriented() const noexcept { return oriented_; }

        //- Return non-const reference to the oriented flag
        orientedOption& oriented() noexcept { return oriented_; }

        //- True if \c ORIENTED
        bool is_oriented() const noexcept { return (oriented_ == ORIENTED); }

        //- Set the oriented flag: on/off
        void setOriented(bool on = true) noexcept
        {
            oriented_ = on ? ORIENTED : UNORIENTED;
        }

        //- Read the "oriented" state from dictionary
        void read(const dictionary& dict);

        //- Write the "oriented" flag entry (if \c ORIENTED)
        //  \return True if entry was written
        bool writeEntry(Ostream& os) const;


    // Member Operators

        void operator+=(const orientedType& ot);
        void operator-=(const orientedType& ot);
        void operator*=(const orientedType& ot);
        void operator/=(const orientedType& ot);
        void operator*=(const scalar s);
        void operator/=(const scalar s);

        //- Convert to bool. Same as is_oriented()
        bool operator()() const noexcept { return is_oriented(); }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, orientedType& ot);
        friend Ostream& operator<<(Ostream& os, const orientedType& ot);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

orientedType min(const orientedType& a, const orientedType& b);
orientedType max(const orientedType& a, const orientedType& b);
orientedType lerp(const orientedType& a, const orientedType& b);

orientedType cmptMultiply(const orientedType& ot1, const orientedType& ot2);
orientedType cmptDivide(const orientedType& ot1, const orientedType& ot);
orientedType cmptAv(const orientedType& ot);

orientedType pow(const orientedType& ot, const scalar p);
orientedType sqr(const orientedType& ot);
orientedType pow2(const orientedType& ot);
orientedType pow3(const orientedType& ot);
orientedType pow4(const orientedType& ot);
orientedType pow5(const orientedType& ot);
orientedType pow6(const orientedType& ot);
orientedType pow025(const orientedType& ot);


orientedType sqrt(const orientedType& ot);
orientedType cbrt(const orientedType& ot);
orientedType magSqr(const orientedType& ot);
orientedType mag(const orientedType& ot);
orientedType sign(const orientedType& ot);
orientedType pos(const orientedType& ot);
orientedType pos0(const orientedType& ot);
orientedType neg(const orientedType& ot);
orientedType neg0(const orientedType& ot);
orientedType posPart(const orientedType& ot);
orientedType negPart(const orientedType& ot);
orientedType inv(const orientedType& ot);


orientedType trans(const orientedType& ot);
orientedType atan2(const orientedType& ot1, const orientedType& ot2);
orientedType hypot(const orientedType& ot1, const orientedType& ot2);
orientedType transform(const orientedType& ot);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

orientedType operator-(const orientedType& ot);
orientedType operator*(const scalar s, const orientedType& ot);
orientedType operator/(const orientedType& ot, const scalar s);

orientedType operator+(const orientedType& ot1, const orientedType& ot2);
orientedType operator-(const orientedType& ot1, const orientedType& ot2);
orientedType operator/(const orientedType& ot1, const orientedType& ot2);
orientedType operator*(const orientedType& ot1, const orientedType& ot2);
orientedType operator^(const orientedType& ot1, const orientedType& ot2);
orientedType operator&(const orientedType& ot1, const orientedType& ot2);
orientedType operator&&(const orientedType& ot1, const orientedType& ot2);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
