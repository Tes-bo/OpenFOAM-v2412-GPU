/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
inline Foam::scalar Foam::Function1Types::Sine<Type>::cycle
(
    const scalar t
) const
{
    // The cycle: (freq * time) or (time / period)
    return
    (
        frequency_
      ? (t - t0_) * frequency_->value(t)
      : (t - t0_) / (period_->value(t) + VSMALL)
    );
}


template<class Type>
inline Foam::scalar
Foam::Function1Types::Sine<Type>::cosForm(const scalar t) const
{
    return
    (
        cos(constant::mathematical::twoPi * this->cycle(t))
      * (amplitude_ ? amplitude_->value(t) : 1.0)
    );
}


template<class Type>
inline Foam::scalar
Foam::Function1Types::Sine<Type>::sinForm(const scalar t) const
{
    return
    (
        sin(constant::mathematical::twoPi * this->cycle(t))
      * (amplitude_ ? amplitude_->value(t) : 1.0)
    );
}


template<class Type>
inline Foam::scalar
Foam::Function1Types::Sine<Type>::squareForm
(
    const scalar t,
    const scalar posFrac
) const
{
    const scalar cyc = this->cycle(t);

    return
    (
        // Fraction of incomplete cycle
        ((cyc - std::floor(cyc)) < posFrac ? 1.0 : -1.0)
      * (amplitude_ ? amplitude_->value(t) : 1.0)
    );
}


template<class Type>
inline Type Foam::Function1Types::Sine<Type>::cosValue(const scalar t) const
{
    return
    (
        cosForm(t) * scale_->value(t) + level_->value(t)
    );
}


template<class Type>
inline Type Foam::Function1Types::Sine<Type>::sinValue(const scalar t) const
{
    return
    (
        sinForm(t) * scale_->value(t) + level_->value(t)
    );
}


template<class Type>
inline Type Foam::Function1Types::Sine<Type>::squareValue
(
    const scalar t,
    const scalar posFrac
) const
{
    return
    (
        squareForm(t, posFrac) * scale_->value(t) + level_->value(t)
    );
}


// ************************************************************************* //
