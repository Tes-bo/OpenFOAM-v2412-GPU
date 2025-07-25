/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

template<class Cmpt>
inline Foam::RowVector<Cmpt>::RowVector(const Foam::zero)
:
    RowVector::msType(Zero)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::RowVector<Cmpt>::RowVector
(
    const MatrixSpace<RowVector<Cmpt2>, Cmpt2, 1, 3>& ms
)
:
    RowVector::msType(ms)
{}


template<class Cmpt>
inline Foam::RowVector<Cmpt>::RowVector
(
    const Cmpt& rvx,
    const Cmpt& rvy,
    const Cmpt& rvz
)
{
    this->v_[X] = rvx;
    this->v_[Y] = rvy;
    this->v_[Z] = rvz;
}


template<class Cmpt>
inline Foam::RowVector<Cmpt>::RowVector(Istream& is)
:
    RowVector::msType(is)
{}


// ************************************************************************* //
