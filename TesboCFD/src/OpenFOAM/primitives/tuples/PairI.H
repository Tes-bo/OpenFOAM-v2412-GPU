/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
inline int Foam::Pair<T>::compare(const Pair<T>& a, const Pair<T>& b)
{
    if (a.first() == b.first() && a.second() == b.second())
    {
        return 1;
    }
    if (a.first() == b.second() && a.second() == b.first())
    {
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Pair<T>::Pair(const T& f, const T& s)
{
    first()  = f;
    second() = s;
}


template<class T>
inline Foam::Pair<T>::Pair(T&& f, T&& s)
{
    first()  = std::move(f);
    second() = std::move(s);
}


template<class T>
inline Foam::Pair<T>::Pair(const std::pair<T,T>& vals)
{
    first()  = vals.first;
    second() = vals.second;
}


template<class T>
inline Foam::Pair<T>::Pair(std::pair<T,T>&& vals)
{
    first()  = std::move(vals.first);
    second() = std::move(vals.second);
}


template<class T>
inline Foam::Pair<T>::Pair(const FixedList<T, 2>& list)
:
    FixedList<T, 2>(list)
{}


template<class T>
inline Foam::Pair<T>::Pair(const T& f, const T& s, const bool doSort)
{
    if (doSort && s < f)
    {
        first()  = s;
        second() = f;
    }
    else
    {
        first()  = f;
        second() = s;
    }
}


template<class T>
inline Foam::Pair<T>::Pair(const FixedList<T, 2>& list, const bool doSort)
:
    Pair<T>(list.template get<0>(), list.template get<1>(), doSort)
{}


template<class T>
inline Foam::Pair<T>::Pair(Istream& is)
:
    FixedList<T, 2>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::Pair<T>::other(const T& val) const
{
    if (first() == second())
    {
        FatalErrorInFunction
            << "Call to other only valid for Pair with differing elements:"
            << *this << abort(FatalError);
    }
    else if (val == first())
    {
        return second();
    }
    else if (val != second())
    {
        FatalErrorInFunction
            << "Pair " << *this
            << " does not contain " << val << abort(FatalError);
    }

    return first();
}


template<class T>
inline void Foam::Pair<T>::flip()
{
    Foam::Swap(first(), second());
}


template<class T>
inline bool Foam::Pair<T>::is_sorted() const
{
    return !(second() < first());
}


template<class T>
inline void Foam::Pair<T>::sort()
{
    if (second() < first())
    {
        flip();
    }
}


// ************************************************************************* //
