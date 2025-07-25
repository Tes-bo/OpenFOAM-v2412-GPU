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
    Foam::predicates::scalars

Description
    A list of unary predicates (tests) on scalars.
    Includes a number of standard comparison predicates
    (eg, "less", "greater", ...)

Note
    This class is still subject to larger changes (2018-11) as it matures.

SeeAlso
    Foam::scalarRange, Foam::scalarRanges

SourceFiles
    scalarPredicates.C
    scalarPredicatesI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarPredicates_H
#define Foam_scalarPredicates_H

#include "scalar.H"
#include "predicates.H"
#include "Enum.H"
#include "List.H"
#include "Tuple2.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class scalars Declaration
\*---------------------------------------------------------------------------*/

class scalars
:
    public List<std::function<bool(Foam::scalar)>>
{
public:

    // Public types

        //- The unary function type for testing a scalar
        typedef std::function<bool(Foam::scalar)> unary;

        //- Enumerations for some standard comparison predicates
        enum opType
        {
            EQUAL,              //!< "eq", "equal"
            NOT_EQUAL,          //!< "neq", "notEqual"
            LESS,               //!< "lt", "less"
            LESS_EQUAL,         //!< "le", "lessEqual", "lessEq"
            GREATER,            //!< "gt", "greater"
            GREATER_EQUAL,      //!< "ge", "greaterEqual", "greaterEq"
            ALWAYS,             //!< Always matches
            NEVER,              //!< Never matches

            // OpenFOAM-v2112 and earlier:
            LESS_EQ = LESS_EQUAL,
            GREATER_EQ = GREATER_EQUAL,
        };

        //- Names for the opType enumeration.
        static const Enum<opType> opNames;


    // Standard comparison methods

        //- Compare for equality, with specified tolerance (non-negative)
        static unary equalOp
        (
            const scalar opVal,
            const scalar tol = VSMALL
        )
        {
            return [=](const scalar val)
            {
                return (Foam::mag(opVal - val) <= tol);
            };
        }


        //- Compare for inequality, with specified tolerance (non-negative)
        static unary notEqualOp
        (
            const scalar opVal,
            const scalar tol = VSMALL
        )
        {
            return [=](const scalar val)
            {
                return (Foam::mag(opVal - val) > tol);
            };
        }


        //- Test if value is 'less' than prescribed
        static unary lessOp(const scalar opVal)
        {
            return std::bind
            (
                std::less<scalar>(), std::placeholders::_1, opVal
            );
        }


        //- Test if value is 'less_equal' to prescribed
        static unary lessEqualOp(const scalar opVal)
        {
            return std::bind
            (
                std::less_equal<scalar>(), std::placeholders::_1, opVal
            );
        }


        //- Test if value is 'greater' than prescribed
        //- Compare scalar values for inequality
        static unary greaterOp(const scalar opVal)
        {
            return std::bind
            (
                std::greater<scalar>(), std::placeholders::_1, opVal
            );
        }


        //- Test if value is 'greater_equal' to prescribed
        static unary greaterEqualOp(const scalar opVal)
        {
            return std::bind
            (
                std::greater_equal<scalar>(), std::placeholders::_1, opVal
            );
        }


    // Special purpose comparison methods

        //- Predicate that always returns true
        static unary trueOp()
        {
            return [](const scalar) { return true; };
        }

        //- Predicate that always returns false
        static unary falseOp()
        {
            return [](const scalar) { return false; };
        }


    // Combining operations

        //- Combine unary tests as an AND operation
        static unary andOp(const unary& test1, const unary& test2)
        {
            return [=](const scalar value)
            {
                return (test1(value) && test2(value));
            };
        }

        //- Combine unary tests as an OR operation
        static unary orOp(const unary& test1, const unary& test2)
        {
            return [=](const scalar value)
            {
                return (test1(value) || test2(value));
            };
        }


    // Factory for standard comparison methods

        //- Standard comparison method by type
        static unary operation
        (
            const opType op,
            const scalar opVal,
            const scalar tol = VSMALL
        );

        //- Standard comparison method by name
        inline static unary operation
        (
            const word& opName,
            const scalar opVal,
            const scalar tol = VSMALL
        );

        //- Standard comparison method by (name, value)
        inline static unary operation
        (
            const Tuple2<word, scalar>& op,
            const scalar tol = VSMALL
        );

        //- Standard comparison method by (name, value)
        inline static unary operation
        (
            const std::pair<word, scalar>& op,
            const scalar tol = VSMALL
        );


    // Constructors

        //- Inherit constructors from List of unary functions
        using List<unary>::List;

        //- Construct from an initializer list of (opName opValue) tuples
        explicit scalars
        (
            std::initializer_list<std::pair<word, scalar>> entries
        );

        //- Copy construct from a list of (opName opValue) tuples
        explicit scalars(const UList<Tuple2<word, scalar>>& entries);

        //- Construct from Istream, from list of (opName opValue) tuples
        explicit scalars(Istream& is);


    //- Destructor
    ~scalars() = default;


    // Member Functions

        //- Assign new predicates from an initializer list of
        //- (opName opValue) tuples
        void assign(std::initializer_list<std::pair<word, scalar>> entries);

        //- Assign new predicates from a list of (opName opValue) tuples
        void assign(const UList<Tuple2<word, scalar>>& entries);


    // Search

        //- Index of the first matching test for the value.
        //  Any occurrences before the start pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label find(const scalar value, label pos = 0) const;

        //- Index of the last match for the value.
        //  Any occurrences after the end pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label rfind(const scalar value, label pos = -1) const;

        //- Matches \em any condition in the list after the offset position.
        //  Linear search from start pos until the end of the list.
        //  Any occurrences before the start pos are ignored.
        //
        //  \return True if the value matches any condition in the (sub)list.
        inline bool contains(const scalar value, label pos = 0) const;

        //- Matches \em any condition in the list.
        //
        //  \return True if the value matches any condition in the list.
        inline bool match(const scalar value) const;

        //- Matches \em any condition in the list.
        //
        //  \return True if the value matches any condition in the list.
        inline bool matchAny(const scalar value) const;

        //- Matches \em all conditions in the list.
        //
        //  \return True if the value matches all conditions in the list.
        inline bool matchAll(const scalar value) const;

        //- Extract list indices for all matches.
        //
        //  \return The locations (indices) of the predicates that match
        inline labelList matching(const scalar value) const;

        //- Extract list indices for all matches.
        //
        //  \param input  A list of scalar inputs to test
        //  \param invert invert the matching logic
        //  \return The locations (indices) in the input list where match()
        //      is true
        inline labelList matching
        (
            const UList<scalar>& input,
            const bool invert=false
        ) const;


    // Member Operators

        //- Identical to contains(), match(), for use as a predicate.
        bool operator()(const scalar value) const { return matchAny(value); }


    // Housekeeping

        //- Same as contains()
        bool found(const scalar value, label pos = 0) const
        {
            return this->contains(value, pos);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates


// IOstream Operators

//- Suppressed write
Ostream& operator<<(Ostream& os, const predicates::scalars& list)
{
    return os;
}

//- Read list of (opName opValue) tuple
Istream& operator>>(Istream& is, predicates::scalars& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarPredicatesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
