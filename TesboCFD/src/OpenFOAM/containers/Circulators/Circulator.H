/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::Circulator

Description
    Walks over a container as if it were circular. The container must have the
    following members defined:
        - size_type
        - difference_type
        - iterator / const_iterator
        - reference / const_reference

    Examples

    \code
        face f(identity(5));

        // Construct Circulator from the face
        Circulator<face> circ(f);

        // If it has a size to iterate over,
        // circulate around the list starting and finishing at the fulcrum.
        if (!circ.empty())
        {
            do
            {
                *circ += 1;

                Info<< "Iterate forwards over face : " << *circ << endl;

            } while (circ.circulate(CirculatorBase::FORWARD));
        }
    \endcode

SourceFiles
    CirculatorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Circulator_H
#define Foam_Circulator_H

#include <type_traits>  // std::conditional

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CirculatorBase Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant parts for circulators
class CirculatorBase
{
public:

    // Public Data

        //- Direction type enumeration
        enum direction
        {
            NONE,
            CLOCKWISE,
            ANTICLOCKWISE
        };


    // Constructors

        //- Default construct
        CirculatorBase() = default;
};


/*---------------------------------------------------------------------------*\
                       Class CirculatorIters Declaration
\*---------------------------------------------------------------------------*/

//- A pair of begin/end iterators used for implementing circular iteration
template<class Container, bool Const>
class CirculatorIters
:
    public CirculatorBase
{
public:

    // STL type definitions

        //- The type that can represent the size of Container
        using size_type = typename Container::size_type;

        //- The type that represents difference between iterator objects
        using difference_type = typename Container::difference_type;

        //- The container iterator type (const/non-const)
        using iterator = typename std::conditional
        <
            Const,
            typename Container::const_iterator,
            typename Container::iterator
        >::type;

        //- The reference type (const/non-const)
        using reference = typename std::conditional
        <
            Const,
            typename Container::const_reference,
            typename Container::reference
        >::type;


private:

    // Private Data

        //- The container beginning
        iterator begin_;

        //- The container end
        iterator end_;

        //- Random access iterator for traversing Container.
        iterator iter_;

        //- Iterator holding location of the fulcrum (start and end) of
        //- the container.
        //- Used to decide when iterator should stop circulating the container
        iterator fulcrum_;


protected:

    // Protected Member Functions

        //- Compare for equality
        inline bool equal(const CirculatorIters<Container, Const>& rhs);

        //- Move iterator forward
        inline void increment();

        //- Move iterator backward
        inline void decrement();


public:

    // Constructors

        //- Default construct
        inline CirculatorIters();

        //- Construct from begin/end iterators
        inline CirculatorIters(const iterator& begin, const iterator& end);

        //- Copy construct
        CirculatorIters(const CirculatorIters<Container, Const>& rhs);


    // Member Functions

        //- True if begin/end iterators are identical
        inline bool empty() const;

        //- Return the range of the iterator pair
        inline size_type size() const;

        //- The distance between the iterator and the fulcrum.
        //  This is equivalent to the number of rotations of the circulator
        inline difference_type nRotations() const;

        //- Circulate around the list in the given direction
        //- \return True if iterator has not yet reached the fulcrum
        inline bool circulate
        (
            const CirculatorBase::direction dir = CirculatorBase::NONE
        );

        //- Set the fulcrum to the current position of the iterator
        inline void setFulcrumToIterator();

        //- Set the iterator to the current position of the fulcrum
        inline void setIteratorToFulcrum();

        //- Dereference the current iterator
        inline reference curr() const;

        //- Dereference the next iterator
        inline reference next() const;

        //- Dereference the previous iterator
        inline reference prev() const;


    // Member Operators

        //- Assignment operator for circulators operating
        //- on the same container type
        inline void operator=(const CirculatorIters<Container, Const>& rhs);

        //- Prefix increment the iterator.
        //  Wraps iterator to the begin if it reaches the end
        inline CirculatorIters<Container, Const>& operator++();

        //- Postfix increment the iterator.
        //  Wraps iterator to the begin if it reaches the end
        inline CirculatorIters<Container, Const> operator++(int);

        //- Prefix decrement the iterator.
        //  Wraps iterator to the end if it reaches the beginning
        inline CirculatorIters<Container, Const>& operator--();

        //- Postfix decrement the iterator.
        //  Wraps iterator to the end if it reaches the beginning
        inline CirculatorIters<Container, Const> operator--(int);

        //- Check for equality of this iterator with another iterator that
        //- operate on the same container type
        inline bool operator==
        (
            const CirculatorIters<Container, Const>&
        ) const;

        //- Check for inequality of this iterator with another iterator that
        //- operate on the same container type
        inline bool operator!=
        (
            const CirculatorIters<Container, Const>&
        ) const;

        //- Dereference the iterator. Same as curr()
        inline reference operator*() const;

        //- Dereference the iterator. Same as curr()
        inline reference operator()() const;

        //- Return the difference between this iterator and another iterator
        //- that operate on the same container type
        inline difference_type operator-
        (
            const CirculatorIters<Container, Const>&
        ) const;
};


/*---------------------------------------------------------------------------*\
                         Class Circulator Declaration
\*---------------------------------------------------------------------------*/

template<class Container>
class Circulator
:
    public CirculatorIters<Container, false>
{
public:

    // Constructors

        //- Default construct
        Circulator() = default;

        //- Construct from begin/end of a container
        explicit Circulator(Container& obj)
        :
            CirculatorIters<Container, false>(obj.begin(), obj.end())
        {}

        //- Construct from two iterators
        Circulator
        (
            const typename Container::iterator& begin,
            const typename Container::iterator& end
        )
        :
            CirculatorIters<Container, false>(begin, end)
        {}

        //- Copy construct
        Circulator(const Circulator<Container>& rhs) = default;


    // Member Operators

        //- Copy assignment
        Circulator<Container>& operator=
        (
            const Circulator<Container>& rhs
        ) = default;
};


/*---------------------------------------------------------------------------*\
                      Class ConstCirculator Declaration
\*---------------------------------------------------------------------------*/

//- Like Foam::Circulator, but with const-access iterators
template<class Container>
class ConstCirculator
:
    public CirculatorIters<Container, true>
{
public:

    // Constructors

        //- Default construct
        ConstCirculator() = default;

        //- Construct from begin/end of a container
        explicit ConstCirculator(const Container& obj)
        :
            CirculatorIters<Container, true>(obj.begin(), obj.end())
        {}

        //- Construct from two iterators
        inline ConstCirculator
        (
            const typename Container::const_iterator& begin,
            const typename Container::const_iterator& end
        )
        :
            CirculatorIters<Container, true>(begin, end)
        {}

        //- Copy construct
        ConstCirculator(const ConstCirculator<Container>& rhs) = default;


    // Member Operators

        //- Copy assignment
        ConstCirculator<Container>& operator=
        (
            const ConstCirculator<Container>& rhs
        ) = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CirculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
