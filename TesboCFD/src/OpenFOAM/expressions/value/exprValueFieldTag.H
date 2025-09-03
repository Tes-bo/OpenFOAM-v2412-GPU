/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::expressions::exprValueFieldTag

Description
    An expressions::exprValue (polymorphic typed union) with an additional
    flag for tracking Field content as uniform etc.

SourceFiles
    exprValueFieldTag.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprValueFieldTag_H
#define Foam_expressions_exprValueFieldTag_H

#include "exprValue.H"
#include "UList.H"  // For ListPolicy

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprValueFieldTag Declaration
\*---------------------------------------------------------------------------*/

class exprValueFieldTag
{
    // Private Data

        //- Uniformity of field (0: empty, 1: uniform, 2: non-uniform, ...)
        //  Values as per internal enum Foam::Detail::ListPolicy::uniformity
        int uniformity_{};

        //- Representative (uniform) type/value for the field
        expressions::exprValue value_{};


public:

    // Generated Methods

        //- Default construct. Uniformity = empty, type = none
        exprValueFieldTag() = default;

        //- Copy construct
        exprValueFieldTag(const exprValueFieldTag&) = default;

        //- Copy assignment
        exprValueFieldTag& operator=(const exprValueFieldTag&) = default;


    // Constructors

        //- Construct as uniform with the specified type/value
        template<class Type>
        explicit exprValueFieldTag(const Type& val)
        {
            set_value(val);
        }

        //- Construct from a range of values
        template<class Type>
        explicit exprValueFieldTag(const Type* first, const Type* last)
        {
            set(first, last);
        }


    // Factory Methods

        //- Make an empty field tag with the specified type (zero-value)
        template<class Type>
        static exprValueFieldTag make_empty()
        {
            exprValueFieldTag tag;  // construct empty, no type
            tag.set_zero<Type>();   // set type and zero value
            return tag;
        }


    // Member Functions

        //- True if the uniformity is "empty"
        bool empty() const noexcept;

        //- True if the uniformity is "uniform"
        bool is_uniform() const noexcept;

        //- True if the uniformity is "non-uniform"
        bool is_nonuniform() const noexcept;

        //- Representative (uniform) value for the field
        const expressions::exprValue& value() const noexcept;

        //- Compare (uniformity, type, value)
        int compare(const exprValueFieldTag& rhs) const;

        //- Test for equality of the values
        bool equal(const exprValueFieldTag& rhs) const;


    // Setters

        //- Set value and uniformity from range of data
        template<class Type>
        void set(const Type* first, const Type* last)
        {
            uniformity_ =
                Foam::Detail::ListPolicy::check_uniformity(first, last);

            if (first != last)  // or uniformity_ != EMPTY
            {
                value_.set<Type>(*first);
            }
            else
            {
                // Is empty, set zero value
                value_.set<Type>(Foam::zero{});
            }
        }

        //- Set as uniform, with specified value
        template<class Type>
        void set_value(const Type& val)
        {
            uniformity_ = Foam::Detail::ListPolicy::uniformity::UNIFORM;
            value_.set<Type>(val);
        }

        //- Set type and zero value (does not affect uniformity)
        template<class Type>
        void set_zero()
        {
            value_.set<Type>(Foam::zero{});
        }

        //- Set as empty with zero value, leave type unchanged
        void set_empty();

        //- Set as non-uniform with zero value, leave type unchanged
        void set_nouniform();


    // Parallel

        //- Inplace parallel reduction, uses worldComm
        void reduce();

        //- Perform a reduction on a copy and return the result
        static exprValueFieldTag returnReduce(const exprValueFieldTag& tag);


    // Reduction operations

        //- Inplace combine - eg, for global uniformity
        void combine(const exprValueFieldTag& b);

        //- Binary combine operator, e.g. for global reduction
        struct combineOp
        {
            exprValueFieldTag operator()
            (
                const exprValueFieldTag& a,
                const exprValueFieldTag& b
            ) const
            {
                exprValueFieldTag result(a);
                result.combine(b);
                return result;
            }
        };


    // IO Operations

        //- Read uniformity label and the value as pair
        void read(Istream& is);

        //- Write uniformity label and the value as pair
        void write(Ostream& os) const;

        //- Print description to Ostream
        void print(Ostream& os) const;


    // Member Operators

        //- Assign from zero. Changes value but not type
        void operator=(const Foam::zero) { value_ = Foam::zero{}; }

        //- Compare (uniformity,value) for equality
        bool operator==(const exprValueFieldTag&) const;

        //- Compare (uniformity,value) for inequality
        bool operator!=(const exprValueFieldTag& rhs) const
        {
            return !(*this == rhs);
        }

        //- Compare (uniformity,value)
        bool operator<(const exprValueFieldTag&) const;
};

} // End namespace expressions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IO
Istream& operator>>(Istream&, expressions::exprValueFieldTag&);
Ostream& operator<<(Ostream&, const expressions::exprValueFieldTag&);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
