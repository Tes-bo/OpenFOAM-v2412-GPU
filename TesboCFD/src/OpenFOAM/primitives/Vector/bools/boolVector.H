/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::boolVector

Description
    Specialized bundling of boolean values as a vector of 3 components,
    element access using x(), y() and z() member functions.
    It also has some methods similar to bitSet.

Note
    The boolVector is not derived from Vector or VectorSpace since
    it does not share very many vector-like characteristics.

SourceFiles
    boolVectorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_boolVector_H
#define Foam_boolVector_H

#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class boolVector Declaration
\*---------------------------------------------------------------------------*/

class boolVector
:
    public FixedList<bool, 3>
{
public:

    // Typedefs

        //- The component type is bool
        typedef bool cmptType;


    // Member Constants

        //- Rank of a vector is 1
        static constexpr direction rank = 1;

        //- Number of components in this vector space
        static constexpr direction nComponents = 3;

        //- Component labeling enumeration
        enum components { X, Y, Z };


    // Generated Methods

        //- Copy construct
        boolVector(const boolVector&) = default;

        //- Copy assignment
        boolVector& operator=(const boolVector&) = default;

        //- Move construct
        boolVector(boolVector&&) = default;

        //- Move assignment
        boolVector& operator=(boolVector&&) = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Default construct, zero-initialized (ie, false)
        inline boolVector();

        //- Uniform construct with specified value
        inline explicit boolVector(const bool val);

        //- Construct from three components
        inline boolVector(const bool vx, const bool vy, const bool vz);

        //- Construct from Istream
        inline explicit boolVector(Istream& is);


    // Member Functions

    // Query

        //- True if all components are set
        //
        //  \note Method name compatibility with bitSet
        inline bool all() const noexcept;

        //- True if any components are set
        //
        //  \note Method name compatibility with bitSet
        inline bool any() const noexcept;

        //- True if no components are set
        //
        //  \note Method name compatibility with bitSet
        inline bool none() const noexcept;

        //- Count number of items set.
        //  \param on can be set to false to count the number of unset bits
        //     instead.
        //
        //  \note Method name compatibility with bitSet
        inline unsigned int count(const bool on=true) const;


    // Component Access

        //- The x component
        bool x() const noexcept { return get<0>(); }

        //- The y component
        bool y() const noexcept { return get<1>(); }

        //- The z component
        bool z() const noexcept { return get<2>(); }

        //- The x component
        bool& x() noexcept { return get<0>(); }

        //- The y component
        bool& y() noexcept { return get<1>(); }

        //- The z component
        bool& z() noexcept { return get<2>(); }


    // Edit

        //- Invert all values
        //
        //  \note Method name compatibility with bitSet
        inline void flip();


    // Operators

        //- Assignment of all entries to the given value
        inline void operator=(const bool value);
};

// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- A boolVector is contiguous (FixedList of bool)
template<> struct is_contiguous<boolVector> : std::true_type {};


/*---------------------------------------------------------------------------*\
                      Specialization pTraits<boolVector>
\*---------------------------------------------------------------------------*/

//- A boolVector has bool data components.
template<>
struct pTraits_cmptType<boolVector> { typedef bool type; };

//- A boolVector has three data components
template<>
struct pTraits_nComponents<boolVector>
:
    std::integral_constant<Foam::direction, 3>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boolVectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
