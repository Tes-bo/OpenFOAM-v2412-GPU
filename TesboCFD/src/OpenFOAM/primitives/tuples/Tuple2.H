/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::Tuple2

Description
    A 2-tuple for storing two objects of dissimilar types.
    The container is similar in purpose to std::pair, but does not expose
    its members directly.

See also
    Foam::Pair for storing two objects of identical types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Tuple2_H
#define Foam_Tuple2_H

#include "Istream.H"
#include "Ostream.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Tuple2 Declaration
\*---------------------------------------------------------------------------*/

template<class T1, class T2 = T1>
class Tuple2
{
    // Private Data

        T1 f_;
        T2 s_;

public:

    // Typedefs (cf. std::pair)

        //- Type of member first, the first template parameter (T1)
        typedef T1 first_type;

        //- Type of member second, the second template parameter (T2)
        typedef T2 second_type;


    // Constructors

        //- Default construct
        Tuple2() = default;

        //- Copy construct from components
        Tuple2(const T1& f, const T2& s)
        :
            f_(f),
            s_(s)
        {}

        //- Move construct from components
        Tuple2(T1&& f, T2&& s)
        :
            f_(std::move(f)),
            s_(std::move(s))
        {}

        //- Copy construct from std::pair
        Tuple2(const std::pair<T1,T2>& vals)
        :
            f_(vals.first),
            s_(vals.second)
        {}

        //- Move construct from std::pair
        Tuple2(std::pair<T1,T2>&& vals)
        :
            f_(std::move(vals.first)),
            s_(std::move(vals.second))
        {}

        //- Construct from Istream
        explicit Tuple2(Istream& is)
        {
            is.readBegin("Tuple2");
            is >> f_ >> s_;
            is.readEnd("Tuple2");
            is.check(FUNCTION_NAME);
        }


    // Member Functions

        //- Access the first element
        const T1& first() const noexcept { return f_; }

        //- Access the first element
        T1& first() noexcept { return f_; }

        //- Access the second element
        const T2& second() const noexcept { return s_; }

        //- Access the second element
        T2& second() noexcept { return s_; }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Hashing for Tuple2 data
template<class T1, class T2>
struct Hash<Tuple2<T1, T2>>
{
    unsigned operator()(const Tuple2<T1, T2>& obj, unsigned seed=0) const
    {
        return Hash<T2>()(obj.second(), Hash<T1>()(obj.first(), seed));
    }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Tuple2
template<class T1, class T2>
inline Tuple2<T2, T1> reverse(const Tuple2<T1,T2>& t)
{
    return Tuple2<T2, T1>(t.second(), t.first());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

template<class T1, class T2>
inline bool operator==(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class T1, class T2>
inline bool operator!=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(a == b);
}


template<class T1, class T2>
inline bool operator<(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return
    (
        a.first() < b.first()
     || (!(b.first() < a.first()) && a.second() < b.second())
    );
}


template<class T1, class T2>
inline bool operator<=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(b < a);
}


template<class T1, class T2>
inline bool operator>(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return (b < a);
}


template<class T1, class T2>
inline bool operator>=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

// Comparing first only

//- Compare tuple-like containers
//  \return reference to the container with the smaller value of first
template<class T1>
struct minFirstOp
{
    const Pair<T1>& operator()(const Pair<T1>& a, const Pair<T1>& b) const
    {
        return (b.first() < a.first()) ? b : a;
    }

    template<class T2>
    const Tuple2<T1,T2>&
    operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return (b.first() < a.first()) ? b : a;
    }
};


//- Assign tuple-like container to use the one with the smaller value of first
template<class T1>
struct minFirstEqOp
{
    void operator()(Pair<T1>& x, const Pair<T1>& y) const
    {
        if (y.first() < x.first()) x = y;
    }

    template<class T2>
    void operator()(Tuple2<T1,T2>& x, const Tuple2<T1,T2>& y) const
    {
        if (y.first() < x.first()) x = y;
    }
};


//- Compare tuple-like containers
//  \return reference to the container with the larger value of first
template<class T1>
struct maxFirstOp
{
    const Pair<T1>& operator()(const Pair<T1>& a, const Pair<T1>& b) const
    {
        return (a.first() < b.first()) ? b : a;
    }

    template<class T2>
    const Tuple2<T1,T2>&
    operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return (a.first() < b.first()) ? b : a;
    }
};


//- Assign tuple-like container to use the one with the larger value of first
template<class T1>
struct maxFirstEqOp
{
    void operator()(Pair<T1>& x, const Pair<T1>& y) const
    {
        if (x.first() < y.first()) x = y;
    }

    template<class T2>
    void operator()(Tuple2<T1,T2>& x, const Tuple2<T1,T2>& y) const
    {
        if (x.first() < y.first()) x = y;
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read Tuple2 from Istream
template<class T1, class T2>
inline Istream& operator>>(Istream& is, Tuple2<T1,T2>& t)
{
    is.readBegin("Tuple2");
    is >> t.first() >> t.second();
    is.readEnd("Tuple2");

    is.check(FUNCTION_NAME);
    return is;
}


//- Write Tuple2 to Ostream
template<class T1, class T2>
inline Ostream& operator<<(Ostream& os, const Tuple2<T1,T2>& t)
{
    os  << token::BEGIN_LIST
        << t.first() << token::SPACE << t.second()
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
