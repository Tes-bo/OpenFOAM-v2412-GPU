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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::BiIndirectList<T>::BiIndirectList
(
    const UList<T>& posList,
    const UList<T>& negList,
    const labelUList& addr
)
:
    posList_(const_cast<UList<T>&>(posList)),
    negList_(const_cast<UList<T>&>(negList)),
    addr_(addr)
{}


template<class T>
inline Foam::BiIndirectList<T>::BiIndirectList
(
    const UList<T>& posList,
    const UList<T>& negList,
    labelList&& addr
)
:
    posList_(const_cast<UList<T>&>(posList)),
    negList_(const_cast<UList<T>&>(negList)),
    addr_(std::move(addr))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::BiIndirectList<T>::resetAddressing(const labelUList& addr)
{
    addr_ = addr;
}


template<class T>
inline void Foam::BiIndirectList<T>::resetAddressing(labelList&& addr)
{
    addr_.transfer(addr);
}


template<class T>
inline Foam::List<T> Foam::BiIndirectList<T>::list() const
{
    List<T> result(size());

    forAll(*this, i)
    {
        result[i] = operator[](i);
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::BiIndirectList<T>::operator[](const label i)
{
    const label index = addr_[i];

    return (index >= 0 ? posList_[index] : negList_[-index-1]);
}


template<class T>
inline const T& Foam::BiIndirectList<T>::operator[](const label i) const
{
    const label index = addr_[i];

    return (index >= 0 ? posList_[index] : negList_[-index-1]);
}


template<class T>
inline void Foam::BiIndirectList<T>::operator=(const UList<T>& ae)
{
    if (addr_.size() != ae.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addr_.size() << " " << ae.size()
            << abort(FatalError);
    }

    forAll(addr_, i)
    {
        operator[](i) = ae[i];
    }
}


template<class T>
inline void Foam::BiIndirectList<T>::operator=(const T& val)
{
    forAll(addr_, i)
    {
        operator[](i) = val;
    }
}


// ************************************************************************* //
