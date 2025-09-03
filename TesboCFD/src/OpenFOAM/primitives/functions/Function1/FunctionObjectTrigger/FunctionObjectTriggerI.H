/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "FunctionObjectTrigger.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::Function1Types::FunctionObjectTrigger<Type>::active() const
{
    if (triggers_.empty())
    {
        return false;
    }

    const label idx = this->time().functionObjects().triggerIndex();

    return
    (
        idx == labelMin ? defaultValue_ : triggers_.found(idx)
    );
}


template<class Type>
inline const Foam::labelList&
Foam::Function1Types::FunctionObjectTrigger<Type>::triggers() const noexcept
{
    return triggers_;
}


template<class Type>
inline void Foam::Function1Types::FunctionObjectTrigger<Type>::resetTriggers
(
    const labelUList& indices
)
{
    triggers_ = indices;
}


template<class Type>
inline Type Foam::Function1Types::FunctionObjectTrigger<Type>::value
(
    const scalar /*unused*/
) const
{
    return this->active() ? pTraits<Type>::one : pTraits<Type>::zero;
}


template<class Type>
inline Type Foam::Function1Types::FunctionObjectTrigger<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return pTraits<Type>::zero;
}


// ************************************************************************* //
