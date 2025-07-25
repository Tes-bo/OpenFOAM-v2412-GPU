/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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

#include "UniformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline UniformField<Type> min
(
    const UniformField<Type>& a,
    const UniformField<Type>& b
)
{
    return UniformField<Type>(min(a.value(), b.value()));
}


template<class Type, class OtherType>
inline OtherType min(const UniformField<Type>& u, const OtherType& o)
{
    return min(u.value(), o);
}


template<class Type, class OtherType>
inline OtherType min(const OtherType& o, const UniformField<Type>& u)
{
    return min(o, u.value());
}


template<class Type>
inline UniformField<Type> max
(
    const UniformField<Type>& a,
    const UniformField<Type>& b
)
{
    return UniformField<Type>(max(a.value(), b.value()));
}


template<class Type, class OtherType>
inline OtherType max(const UniformField<Type>& u, const OtherType& o)
{
    return max(u.value(), o);
}


template<class Type, class OtherType>
inline OtherType max(const OtherType& o, const UniformField<Type>& u)
{
    return max(o, u.value());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
