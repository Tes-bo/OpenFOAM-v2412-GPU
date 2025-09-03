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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, class Addr>
template<class ListType>
inline void Foam::IndirectListBase<T, Addr>::copyList(const ListType& rhs)
{
    if
    (
        this
     == reinterpret_cast<IndirectListBase<T,Addr>*>(const_cast<ListType*>(&rhs))
    )
    {
        return;  // Self-assignment is a no-op
    }

    const label len = addr_.size();

    if (len != rhs.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: " << len << " " << rhs.size()
            << abort(FatalError);
    }

    // Or std::copy(rhs.cbegin(), rhs.cend(), this->begin());
    for (label i = 0; i < len; ++i)
    {
        values_[addr_[i]] = rhs[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Addr>
inline Foam::IndirectListBase<T, Addr>::IndirectListBase
(
    const UList<T>& values,
    const Addr& addr
)
:
    values_(const_cast<UList<T>&>(values)),
    addr_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Addr>
inline bool Foam::IndirectListBase<T, Addr>::uniform() const
{
    const label len = this->size();

    if (!len)
    {
        return false;
    }

    const T& val = (*this)[0];  // first

    for (label i = 1; i < len; ++i)
    {
        if (val != (*this)[i])
        {
            return false;
        }
    }

    return true;
}


template<class T, class Addr>
inline bool Foam::IndirectListBase<T, Addr>::contains
(
    const T& val,
    label pos,
    label len
) const
{
    return (this->find(val, pos, len) >= 0);
}


template<class T, class Addr>
inline Foam::label Foam::IndirectListBase<T, Addr>::fcIndex(const label i) const
{
    return (i == addr_.size()-1 ? 0 : i+1);
}


template<class T, class Addr>
inline Foam::label Foam::IndirectListBase<T, Addr>::rcIndex(const label i) const
{
    return (i ? i-1 : addr_.size()-1);
}


template<class T, class Addr>
inline const T& Foam::IndirectListBase<T, Addr>::front() const
{
    return values_[addr_.front()];
}

template<class T, class Addr>
inline T& Foam::IndirectListBase<T, Addr>::front()
{
    return values_[addr_.front()];
}


template<class T, class Addr>
inline const T& Foam::IndirectListBase<T, Addr>::back() const
{
    return values_[addr_.back()];
}

template<class T, class Addr>
inline T& Foam::IndirectListBase<T, Addr>::back()
{
    return values_[addr_.back()];
}


template<class T, class Addr>
inline const T& Foam::IndirectListBase<T, Addr>::fcValue(const label i) const
{
    return values_[this->fcIndex(i)];
}


template<class T, class Addr>
inline T& Foam::IndirectListBase<T, Addr>::fcValue(const label i)
{
    return values_[this->fcIndex(i)];
}


template<class T, class Addr>
inline const T& Foam::IndirectListBase<T, Addr>::rcValue(const label i) const
{
    return values_[this->rcIndex(i)];
}


template<class T, class Addr>
inline T& Foam::IndirectListBase<T, Addr>::rcValue(const label i)
{
    return values_[this->rcIndex(i)];
}

template<class T, class Addr>
inline Foam::List<T> Foam::IndirectListBase<T, Addr>::list() const
{
    const label len = addr_.size();

    List<T> result(len);

    // Or std::copy(this->cbegin(), this->cend(), result.begin());
    for (label i = 0; i < len; ++i)
    {
        result[i] = values_[addr_[i]];
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Addr>
inline T& Foam::IndirectListBase<T, Addr>::operator[](const label i)
{
    return values_[addr_[i]];
}


template<class T, class Addr>
inline const T&
Foam::IndirectListBase<T, Addr>::operator[](const label i) const
{
    return values_[addr_[i]];
}


template<class T, class Addr>
inline void Foam::IndirectListBase<T, Addr>::operator=(const T& val)
{
    // Or std::fill(this->begin(), this->end(), val);
    for (const label idx : addr_)
    {
        values_[idx] = val;
    }
}


template<class T, class Addr>
inline void Foam::IndirectListBase<T, Addr>::operator=(const Foam::zero)
{
    // Or std::fill(this->begin(), this->end(), Zero);
    for (const label idx : addr_)
    {
        values_[idx] = Zero;
    }
}


template<class T, class Addr>
inline void Foam::IndirectListBase<T, Addr>::operator=
(
    const UList<T>& rhs
)
{
    this->copyList(rhs);
}


template<class T, class Addr>
inline void Foam::IndirectListBase<T, Addr>::operator=
(
    const IndirectListBase<T, Addr>& rhs
)
{
    this->copyList(rhs);
}


template<class T, class Addr>
template<class AnyAddr>
inline void Foam::IndirectListBase<T, Addr>::operator=
(
    const IndirectListBase<T, AnyAddr>& rhs
)
{
    this->copyList(rhs);
}


// ************************************************************************* //
