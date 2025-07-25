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
    Foam::Pair

Description
    An ordered pair of two objects of type \<T\> with first() and second()
    elements.

SourceFiles
    PairI.H

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Pair_H
#define Foam_Pair_H

#include "FixedList.H"
#include "Istream.H"
#include "Ostream.H"
#include <utility>  // For std::move, std::pair

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class Pair;

// Common pair types
typedef Pair<label> labelPair;      //!< A pair of labels
typedef Pair<word> wordPair;        //!< A pair of words


/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Pair
:
    public FixedList<T, 2>
{
public:

    // Generated Methods

        //- Default construct
        Pair() = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Copy construct from components
        inline Pair(const T& f, const T& s);

        //- Move construct from components
        inline Pair(T&& f, T&& s);

        //- Copy construct from std::pair
        inline Pair(const std::pair<T,T>& vals);

        //- Move construct from std::pair
        inline Pair(std::pair<T,T>&& vals);

        //- Copy construct FixedList of two items
        inline Pair(const FixedList<T, 2>& list);

        //- Copy construct, optionally sorted with first less-than second
        inline Pair(const T& f, const T& s, const bool doSort);

        //- Copy construct, optionally sorted with first less-than second
        inline Pair(const FixedList<T, 2>& list, const bool doSort);

        //- Construct from Istream
        inline explicit Pair(Istream& is);


    // Member Functions

    // Access

        //- Access the first element
        const T& first() const noexcept { return this->template get<0>(); }

        //- Access the first element
        T& first() noexcept { return this->template get<0>(); }

        //- Access the second element
        const T& second() const noexcept { return this->template get<1>(); }

        //- Access the second element
        T& second() noexcept { return this->template get<1>(); }

        //- Return other element
        inline const T& other(const T& a) const;


    // Queries

        //- True if first() is less-than-equal second()
        inline bool is_sorted() const;


    // Editing

        //- Flip the Pair in-place.
        inline void flip();

        //- Sort so that first() is less-than second()
        inline void sort();


    // Comparison

        //- Compare Pairs
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in reversed order
        static inline int compare(const Pair<T>& a, const Pair<T>& b);


    // Hashing

        //- Symmetric hashing functor for Pair, hashes lower value first
        //  Regular hasher inherited from FixedList
        struct symmHasher
        {
            unsigned operator()(const Pair<T>& obj, unsigned seed=0) const
            {
                Foam::Hash<T> op;
                if (obj.second() < obj.first())
                {
                    return op(obj.first(), op(obj.second(), seed));
                }
                else
                {
                    return op(obj.second(), op(obj.first(), seed));
                }
            }
        };


    // Housekeeping

        //- Deprecated(2023-07) Use is_sorted() method
        //  \deprecated(2023-07) Use is_sorted() method
        bool sorted() const { return is_sorted(); }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Pair is contiguous if the type is contiguous
template<class T>
struct is_contiguous<Pair<T>> : is_contiguous<T> {};

//- Check for Pair of labels
template<class T>
struct is_contiguous_label<Pair<T>> : is_contiguous_label<T> {};

//- Check for Pair of scalars
template<class T>
struct is_contiguous_scalar<Pair<T>> : is_contiguous_scalar<T> {};

//- Hashing for Pair of data
template<class T>
struct Hash<Pair<T>> : Pair<T>::hasher {};

//- Hashing for std::pair data
template<class T1, class T2>
struct Hash<std::pair<T1, T2>>
{
    unsigned operator()(const std::pair<T1, T2>& obj, unsigned seed=0) const
    {
        return Hash<T2>()(obj.second, Hash<T1>()(obj.first, seed));
    }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Pair
template<class T>
Pair<T> reverse(const Pair<T>& p)
{
    return Pair<T>(p.second(), p.first());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

template<class T>
bool operator==(const Pair<T>& a, const Pair<T>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class T>
bool operator!=(const Pair<T>& a, const Pair<T>& b)
{
    return !(a == b);
}


template<class T>
bool operator<(const Pair<T>& a, const Pair<T>& b)
{
    return
    (
        a.first() < b.first()
     || (!(b.first() < a.first()) && a.second() < b.second())
    );
}


template<class T>
bool operator<=(const Pair<T>& a, const Pair<T>& b)
{
    return !(b < a);
}


template<class T>
bool operator>(const Pair<T>& a, const Pair<T>& b)
{
    return (b < a);
}


template<class T>
bool operator>=(const Pair<T>& a, const Pair<T>& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read std::pair from Istream
template<class T1, class T2>
inline Istream& operator>>(Istream& is, std::pair<T1,T2>& t)
{
    is.readBegin("pair");
    is >> t.first >> t.second;
    is.readEnd("pair");

    is.check(FUNCTION_NAME);
    return is;
}


//- Write std::pair to Ostream
template<class T1, class T2>
inline Ostream& operator<<(Ostream& os, const std::pair<T1,T2>& t)
{
    os  << token::BEGIN_LIST
        << t.first << token::SPACE << t.second
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
