/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANB WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D(const Foam::zero)
:
    Barycentric2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc
)
{
    this->v_[A] = va;
    this->v_[B] = vb;
    this->v_[C] = vc;
}


template<class Cmpt>
inline Foam::Barycentric2D<Cmpt>::Barycentric2D
(
    const Cmpt& va,
    const Cmpt& vb
)
:
    Barycentric2D<Cmpt>(va, vb, (pTraits<Cmpt>::one - va - vb))
{}


// * * * * * * * * * * * * * * * Member Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt Foam::Barycentric2D<Cmpt>::inner
(
    const Barycentric2D<Cmpt>& b2
) const
{
    const Barycentric2D<Cmpt>& b1 = *this;

    return (b1.a()*b2.a() + b1.b()*b2.b() + b1.c()*b2.c());
}


template<class Cmpt>
inline bool Foam::Barycentric2D<Cmpt>::outside() const
{
    return ((this->v_[A] < 0) || (this->v_[B]) < 0 || (this->v_[C] < 0));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt operator&
(
    const Barycentric2D<Cmpt>& b1,
    const Barycentric2D<Cmpt>& b2
)
{
    return b1.inner(b2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
