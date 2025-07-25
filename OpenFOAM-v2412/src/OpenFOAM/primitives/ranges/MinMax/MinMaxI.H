/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
inline Foam::MinMax<T> Foam::MinMax<T>::ge(const T& minVal)
{
    return MinMax<T>(minVal, pTraits<T>::max);
}


template<class T>
inline Foam::MinMax<T> Foam::MinMax<T>::le(const T& maxVal)
{
    return MinMax<T>(pTraits<T>::min, maxVal);
}


template<class T>
inline Foam::MinMax<T> Foam::MinMax<T>::zero_one()
{
    return MinMax<T>(Foam::zero_one{});
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::MinMax<T>::MinMax()
:
    Tuple2<T,T>(pTraits<T>::max, pTraits<T>::min)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const T& minVal, const T& maxVal)
:
    Tuple2<T,T>(minVal, maxVal)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const std::pair<T,T>& range)
:
    Tuple2<T,T>(range.first, range.second)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const Pair<T>& range)
:
    Tuple2<T,T>(range.first(), range.second())
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const Foam::zero)
:
    Tuple2<T,T>(pTraits<T>::zero, pTraits<T>::zero)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const Foam::zero_one)
:
    Tuple2<T,T>(pTraits<T>::zero, pTraits<T>::one)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const T& val)
:
    Tuple2<T,T>(val, val)
{}


template<class T>
inline Foam::MinMax<T>::MinMax(const UList<T>& vals)
:
    MinMax<T>()
{
    add(vals);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::MinMax<T>::min() const noexcept
{
    return this->first();
}


template<class T>
inline T& Foam::MinMax<T>::min() noexcept
{
    return this->first();
}


template<class T>
inline const T& Foam::MinMax<T>::max() const noexcept
{
    return this->second();
}


template<class T>
inline T& Foam::MinMax<T>::max() noexcept
{
    return this->second();
}


template<class T>
inline T Foam::MinMax<T>::centre() const
{
    // Multiply before adding to avoid possible overflow
    return (0.5 * min()) + (0.5 * max());
}


template<class T>
inline T Foam::MinMax<T>::span() const
{
    return (empty() ? Zero : (max() - min()));
}


template<class T>
inline Foam::scalar Foam::MinMax<T>::mag() const
{
    return ::Foam::mag(span());
}


template<class T>
inline bool Foam::MinMax<T>::empty() const
{
    // Is empty/invalid if (max < min), ie, !(min <= max)
    return (max() < min());
}


template<class T>
inline bool Foam::MinMax<T>::good() const
{
    return !empty();
}


template<class T>
inline void Foam::MinMax<T>::reset()
{
    min() = pTraits<T>::max;
    max() = pTraits<T>::min;
}


template<class T>
inline void Foam::MinMax<T>::reset(const T& val)
{
    min() = val;
    max() = val;
}


template<class T>
inline void Foam::MinMax<T>::reset(const T& minVal, const T& maxVal)
{
    min() = minVal;
    max() = maxVal;
}


template<class T>
inline bool Foam::MinMax<T>::intersects(const MinMax<T>& b) const
{
    const MinMax<T>& a = *this;
    return (a.min() < b.max() && b.min() < a.max());
}


template<class T>
inline bool Foam::MinMax<T>::overlaps(const MinMax<T>& b) const
{
    const MinMax<T>& a = *this;
    return (a.min() <= b.max() && b.min() <= a.max());
}


template<class T>
inline int Foam::MinMax<T>::compare(const T& val) const
{
    if (good())
    {
        if (max() < val)
        {
            return -1;
        }
        else if (val < min())
        {
            return 1;
        }
    }

    return 0;
}


template<class T>
inline bool Foam::MinMax<T>::contains(const T& val) const
{
    return good() && !(val < min() || max() < val);
}


template<class T>
inline T Foam::MinMax<T>::clamp(const T& val) const
{
    if (good())
    {
        return Foam::min(Foam::max(val, min()), max());
    }

    return val; // Pass-through
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::add(const MinMax& other)
{
    min() = Foam::min(min(), other.min());
    max() = Foam::max(max(), other.max());
    return *this;
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::add(const T& val)
{
    min() = Foam::min(min(), val);
    max() = Foam::max(max(), val);
    return *this;
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::add(const UList<T>& vals)
{
    for (const T& val : vals)
    {
        add(val);
    }
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::MinMax<T>::operator()(const T& val) const
{
    return contains(val);
}


// Perhaps not entirely useful
//
// template<class T>
// inline Foam::MinMax<T>& Foam::MinMax<T>::operator-=
// (
//     const MinMax<T>& b
// )
// {
//     MinMax<T>& a = *this;
//
//     // Remove overlapping portion, but cannot create a 'hole' in the middle
//
//     if (a.empty() || b.empty())
//     {
//         // Invalid range(s): no-op
//     }
//     else if (a.min() < b.max() && b.min() <= a.min())
//     {
//         // Overlap on the left
//         min() = ::Foam::max(a.min(), b.max());
//     }
//     else if (b.min() < a.max() && a.max() <= b.max())
//     {
//         // Overlap on the right
//         max() = ::Foam::min(a.max(), b.min());
//     }
//
//     return *this;
// }


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator&=(const MinMax<T>& b)
{
    min() = ::Foam::max(min(), b.min());
    max() = ::Foam::min(max(), b.max());

    return *this;
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator+=(const MinMax<T>& b)
{
    return add(b);
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator+=(const T& val)
{
    return add(val);
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator+=(const UList<T>& vals)
{
    return add(vals);
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator*=(const scalar& s)
{
    min() *= s;
    max() *= s;
    return *this;
}


template<class T>
inline Foam::MinMax<T>& Foam::MinMax<T>::operator/=(const scalar& s)
{
    min() /= s;
    max() /= s;
    return *this;
}


// ************************************************************************* //
