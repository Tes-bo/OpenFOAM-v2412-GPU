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

\*---------------------------------------------------------------------------*/

#include "FixedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list
) noexcept
:
    UList<T>(const_cast<T*>(list.cdata()), list.size())
{}


template<class T>
template<unsigned N>
inline Foam::SubList<T>::SubList
(
    const FixedList<T, N>& list
) noexcept
:
    UList<T>(const_cast<T*>(list.cdata()), list.size())
{}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label len
)
:
    UList<T>(const_cast<T*>(list.cdata()), len)
{
    #ifdef FULLDEBUG
    list.checkSize(len);
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label len,
    const label start
)
:
    UList<T>(const_cast<T*>(list.cdata() + start), len)
{
    #ifdef FULLDEBUG
    list.checkRange(start, len);
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const labelRange& range
)
:
    SubList<T>(range.subset0(list.size()), list)
{
    #ifdef FULLDEBUG
    // The subset0() above always produces valid ranges but want to check
    // that the input itself was valid
    list.checkRange(range.start(), range.size());
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const labelRange& range,
    const UList<T>& list
)
:
    SubList<T>(list, range.size(), range.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset(std::nullptr_t) noexcept
{
    UList<T>::shallowCopy(nullptr);
    return *this;
}


template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset
(
    const UList<T>& list
) noexcept
{
    UList<T>::shallowCopy(list);
    return *this;
}


template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset
(
    const UList<T>& list,
    const label len
)
{
    #ifdef FULLDEBUG
    list.checkSize(len);
    #endif

    UList<T>::shallowCopy(const_cast<T*>(list.cdata()), len);
    return *this;
}


template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset
(
    const UList<T>& list,
    const label len,
    const label start
)
{
    #ifdef FULLDEBUG
    list.checkRange(start, len);
    #endif

    UList<T>::shallowCopy
    (
        const_cast<T*>(list.cdata() + start),
        len
    );
    return *this;
}


template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset
(
    const UList<T>& list,
    const labelRange& range
)
{
    #ifdef FULLDEBUG
    // subset0() always produces valid ranges but want to check
    // that the input itself was valid
    list.checkRange(range.start(), range.size());
    #endif

    labelRange clamped(range.subset0(list.size()));

    UList<T>::shallowCopy
    (
        const_cast<T*>(list.cdata() + clamped.start()),
        clamped.size()
    );
    return *this;
}


template<class T>
inline Foam::UList<T>& Foam::SubList<T>::reset
(
    const labelRange& range,
    const UList<T>& list
)
{
    #ifdef FULLDEBUG
    list.checkRange(range.start(), range.size());
    #endif

    UList<T>::shallowCopy
    (
        const_cast<T*>(list.cdata() + range.start()),
        range.size()
    );
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::operator const Foam::List<T>&() const
{
    return *reinterpret_cast<const List<T>*>(this);
}


template<class T>
inline void Foam::SubList<T>::operator=(const SubList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const UList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
template<class Addr>
inline void Foam::SubList<T>::operator=(const IndirectListBase<T, Addr>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::SubList<T>::operator=(const Foam::zero)
{
    UList<T>::operator=(Foam::zero{});
}


// ************************************************************************* //
