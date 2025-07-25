/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::labelRange

Description
    A range or interval of labels defined by a start and a size.

SourceFiles
    labelRange.C
    labelRangeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_labelRange_H
#define Foam_labelRange_H

#include "IntRange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class MinMax;
template<class T> class Pair;

/*---------------------------------------------------------------------------*\
                         Class labelRange Declaration
\*---------------------------------------------------------------------------*/

class labelRange
:
    public IntRange<label>
{
public:

    // STL type definitions

        //- Input iterator with const access
        using IntRange<label>::const_iterator;

        //- Reverse input iterator with const access
        using IntRange<label>::const_reverse_iterator;


    // Static Data Members

        //- Debugging
        static int debug;


    // Generated Methods: copy/move construct, copy/move assignment


    // Constructors

        //- Default construct an empty range (0,0)
        inline constexpr labelRange() noexcept;

        //- Construct a range with specified length, starting at zero (0,len)
        inline explicit constexpr labelRange(const label len) noexcept;

        //- Construct a range from start/length, no checks
        inline labelRange(const label beg, const label len) noexcept;

        //- Construct a range from start/size, enforces non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
        //  \see adjust()
        inline labelRange
        (
            const label beg,
            const label len,
            const bool adjustStart
        ) noexcept;

        //- Construct from a min/max range (both inclusive),
        //- enforces non-negative size.
        //  Passing an invalid range results in an empty labelRange
        explicit labelRange(const MinMax<label>& range) noexcept;

        //- Construct from start (inclusive) and end (exclusive) values,
        //- enforces non-negative size.
        //  Passing an invalid range results in an empty labelRange
        explicit labelRange(const Pair<label>& start_end) noexcept;

        //- Construct from Istream.
        explicit labelRange(Istream& is);


    // Member Functions

    // Edit

        //- Adjust the start to avoid negative indices.
        //  The size is decreased accordingly, but will never become negative.
        //  Eg, adjusting (-10, 15) becomes (0,5).
        //  adjusting (-20, 15) becomes (0,0)
        void adjust() noexcept;

        //- Reset start and length, no checks
        using IntRange<label>::reset;

        //- Reset start and length, enforces non-negative size.
        //  Optionally adjust the start to avoid any negative indices.
        inline void reset
        (
            const label beg,
            const label len,
            const bool adjustStart
        ) noexcept;


    // Other

        //- Return list of labels corresponding to the range.
        //  Same as Foam::identity()
        List<label> labels() const;

        //- Return true if the ranges overlap.
        //  Optional test for ranges that also just touch each other
        bool overlaps(const labelRange& range, bool touches=false) const;

        //- Return a joined range, squashing any gaps in between
        //  A prior overlaps() check can be used to avoid squashing gaps.
        labelRange join(const labelRange& range) const;

        //- Calculate the intersection of the range with another.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset(const labelRange& range) const;

        //- Calculate the intersection with the given start/size range.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset(const label start, const label size) const;

        //- Calculate the intersection with the given 0/size range.
        //  If there is no intersection, it returns an empty range with zero
        //  for start/size.
        labelRange subset0(const label size) const;


    // Housekeeping

        //- Deprecated(2020-09) True if range has size greater than zero
        //
        //  \deprecated(2020-09) - use bool operator - or good()
        bool valid() const noexcept { return good(); }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Conversion/extraction to labelRange operation (functor).
//  Specializations shall provide a corresponding \c operator().
//  For example,
//  \code
//  template<>
//  struct labelRangeOp<polyPatch>
//  {
//      labelRange operator()(const polyPatch& pp) const
//      {
//          return labelRange(pp.start(), pp.size());
//      }
//  };
//  \endcode
template<class> struct labelRangeOp;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangeI.H"

#endif

// ************************************************************************* //
