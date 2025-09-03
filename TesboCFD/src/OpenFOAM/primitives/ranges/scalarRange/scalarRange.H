/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::scalarRange

Description
    Scalar bounds to be used as a unary predicate.

    The bound can be specified as an "MIN:MAX" range, as a "MIN:" or ":MAX"
    bound or simply as a single "VALUE".

    When defined via the parse() method, the special string "none" can be
    used to define an empty (inverse) range.

SeeAlso
    Foam::MinMax
    Foam::predicates::scalars

SourceFiles
    scalarRange.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarRange_H
#define Foam_scalarRange_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class scalarRange;
template<class T> class MinMax;

Ostream& operator<<(Ostream& os, const scalarRange& range);


/*---------------------------------------------------------------------------*\
                         Class scalarRange Declaration
\*---------------------------------------------------------------------------*/

class scalarRange
{
    //- The type of range
    enum rangeTypes : unsigned char
    {
        NONE = 0,   //!< Empty or invalid (inverse) range. Never matches.
        EQ,         //!< Single value equality
        GE,         //!< Greater-than-equal.
        GT,         //!< Greater-than.
        LE,         //!< Less-than-equal.
        LT,         //!< Less-than.
        GE_LE,      //!< An inclusive range.
        ALWAYS      //!< Placeholder range that always matches
    };


    // Private Member Data

        //- The min value of the range
        scalar min_;

        //- The max value of the range
        scalar max_;

        //- The range type
        enum rangeTypes type_;

        //- Construct from components, no checks
        inline constexpr scalarRange
        (
            const rangeTypes type,
            const scalar minVal,
            const scalar maxVal
        ) noexcept;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef scalar value_type;


    // Static Data Members

        //- A range that always matches
        static const scalarRange always;


    // Constructors

        //- Construct an empty (inverse, NONE) range - never matches
        inline constexpr scalarRange() noexcept;

        //- Construct an exact value matcher
        inline explicit constexpr scalarRange(const scalar val) noexcept;

        //- Construct a range from min-value to max-value
        //  Check validity of the range and sets to NONE or EQ if required.
        inline scalarRange(const scalar minVal, const scalar maxVal) noexcept;

        //- Copy construct from a min/max range.
        //  Automatically decides if this is a GE_LE or NONE range
        explicit scalarRange(const MinMax<label>& range) noexcept;

        //- Copy construct from a min/max range.
        //  Automatically decides if this is a GE_LE or NONE range
        explicit scalarRange(const MinMax<scalar>& range) noexcept;


    // Static Constructors

        //- Construct by parsing string content.
        //  A colon (:) is used as a range marker or when specifying
        //  greater-than or less-than bounds.
        //
        //  \note The special string "none" can be used define an empty
        //      (inverse) range
        //
        //  \return True if no parse problems were encountered.
        static bool parse(const std::string& str, scalarRange& range);

        //- Construct by parsing string content.
        //  Any parse problems are emitted as information and the returned
        //  range is of type empty().
        //  \return The parsed range, which is empty() on any problems
        static scalarRange parse(const std::string& str);


        //- A greater-equals bound
        inline static constexpr scalarRange ge(const scalar minVal) noexcept;

        //- A greater-than bound
        inline static constexpr scalarRange gt(const scalar minVal) noexcept;

        //- A greater-equals zero bound
        inline static constexpr scalarRange ge0() noexcept;

        //- A greater-than zero bound
        inline static constexpr scalarRange gt0() noexcept;

        //- A less-equals bound
        inline static constexpr scalarRange le(const scalar maxVal) noexcept;

        //- A less-than bound
        inline static constexpr scalarRange lt(const scalar maxVal) noexcept;

        //- A greater-equals 0, less-equals 1 bound
        inline static constexpr scalarRange zero_one() noexcept;


    // Member Functions

        //- True if range is empty (eg, inverted, NONE)
        inline bool empty() const noexcept;

        //- True if range is non-empty.
        inline bool good() const noexcept;

        //- True if the range bounds represent a single value.
        inline bool single() const noexcept;

        //- The min value of the range.
        scalar min() const noexcept { return min_; }

        //- The max value of the range.
        scalar max() const noexcept { return max_; }

        //- Reset to an empty (inverse, NONE) range.
        inline void reset() noexcept;

        //- Same as reset() - reset to an empty (inverse, NONE) range.
        void clear() noexcept { reset(); }


        //- A representative (average) value for the range.
        //  For GE, LE bounds it is the min/max value, respectively.
        inline scalar value() const;

        //- True if the value is matched by the condition
        inline bool contains(const scalar val) const;

        //- True if the value is matched by the condition
        bool match(const scalar value) const { return contains(value); }


    // Member Operators

        //- For use as a predicate, same as contains(), match()
        bool operator()(const scalar value) const { return contains(value); }

        inline constexpr bool operator==(const scalarRange& rhs) const noexcept;
        inline constexpr bool operator!=(const scalarRange& rhs) const noexcept;


    // IOstream Operators

        //- Print information about the range
        void print(Ostream& os) const;

        //- Print information about the range
        friend Ostream& operator<<(Ostream& os, const scalarRange& range);


    // Housekeeping

        //- Same as good() or !empty()
        bool valid() const noexcept { return good(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarRangeI.H"

#endif

// ************************************************************************* //
