/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

inline Foam::cell::cell(const label sz)
:
    labelList(sz, -1)
{}


inline Foam::cell::cell(const labelUList& list)
:
    labelList(list)
{}


template<unsigned N>
inline Foam::cell::cell(const FixedList<label, N>& list)
:
    labelList(list)
{}


inline Foam::cell::cell(labelList&& list)
:
    labelList(std::move(list))
{}


inline Foam::cell::cell(const labelUList& list, const labelUList& indices)
:
    labelList(list, indices)
{}


template<unsigned N>
inline Foam::cell::cell
(
    const labelUList& list,
    const FixedList<label, N>& indices
)
:
    labelList(list, indices)
{}


inline Foam::cell::cell(Istream& is)
:
    labelList(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::cell::nFaces() const noexcept
{
    return size();
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool Foam::operator!=(const cell& a, const cell& b)
{
    return !(a == b);
}


// ************************************************************************* //
