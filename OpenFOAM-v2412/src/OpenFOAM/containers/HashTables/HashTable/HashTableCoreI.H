/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * const iterator pair * * * * * * * * * * * * * * //

template<class IteratorType, class TableType>
inline Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>
::const_iterator_pair()
:
    size_(0),
    iter_()
{}


template<class IteratorType, class TableType>
inline Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>
::const_iterator_pair
(
    const TableType& tbl
)
:
    size_(tbl.size()),
    iter_(tbl.begin())
{}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::begin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cbegin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::end() const
{
    return IteratorType();
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cend() const
{
    return IteratorType();
}


// ************************************************************************* //
