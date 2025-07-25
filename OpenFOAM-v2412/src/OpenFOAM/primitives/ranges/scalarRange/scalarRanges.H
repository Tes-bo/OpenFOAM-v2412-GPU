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
    Foam::scalarRanges

Description
    A collection of scalar bounds to be used as a unary predicate.

SeeAlso
    Foam::predicates::scalars

SourceFiles
    scalarRanges.C
    scalarRangesI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarRanges_H
#define Foam_scalarRanges_H

#include "scalarRange.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarRanges Declaration
\*---------------------------------------------------------------------------*/

class scalarRanges
:
    public List<scalarRange>
{
public:

    // Constructors

        //- Inherit constructors from List of scalarRange
        using List<scalarRange>::List;

        //- Default construct
        constexpr scalarRanges() noexcept = default;

        //- Construct by parsing string using the parse() static method
        //- with optional reporting if any range fails to parse
        explicit scalarRanges(const std::string& str, bool report = true)
        :
            List<scalarRange>(scalarRanges::parse(str, report))
        {}


    // Static Constructors

        //- Construct by parsing string for scalar ranges
        //- with optional reporting if any range fails to parse.
        //  The individual items are space, comma or semicolon delimited.
        static scalarRanges parse(const std::string& str, bool report = true);


    // Member Functions

        //- True if the value is contained by any of the sub-ranges
        bool contains(const scalar value) const
        {
            for (const scalarRange& range : *this)
            {
                if (range.contains(value))
                {
                    return true;
                }
            }
            return false;
        }

        //- True if the value is matched by any of the sub-ranges
        bool match(const scalar value) const { return contains(value); }


    // Member Operators

        //- For use as a predicate, same as contains()
        bool operator()(const scalar value) const { return contains(value); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
