/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "UniformValueField.H"
#include "Field.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::UniformValueField<Type>::value
(
    const scalar x
) const
{
    auto tfld =
        tmp<Field<Type>>::New
        (
            this->size(),
            uniformValuePtr_->value(x)
        );

    return this->transform(tfld);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::UniformValueField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    auto tfld =
        tmp<Field<Type>>::New
        (
            this->size(),
            uniformValuePtr_->integrate(x1, x2)
        );

    return this->transform(tfld);
}


// ************************************************************************* //
