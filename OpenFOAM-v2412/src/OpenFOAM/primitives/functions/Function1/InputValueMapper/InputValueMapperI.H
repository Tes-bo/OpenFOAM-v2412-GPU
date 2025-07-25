/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "InputValueMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1Types::InputValueMapper<Type>::value
(
    const scalar t
) const
{
    switch (mappingMode_)
    {
        case mappingMode::NONE:
        {
            return value_->value(t);;
        }
        case mappingMode::FUNCTION1:
        {
            return value_->value(mappingValuePtr_->value(t));
        }
        case mappingMode::MINMAX:
        {
            scalar tlim = clamp(t, min_, max_);

            return value_->value(tlim);
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << mappingModeNames_[mappingMode_]
                << ".  Available options are: " << mappingModeNames_.sortedToc()
                << abort(FatalError);
        }
    }

    return pTraits<Type>::zero;
}


template<class Type>
Type Foam::Function1Types::InputValueMapper<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    switch (mappingMode_)
    {
        case mappingMode::NONE:
        {
            return value_->integrate(x1, x2);
        }
        case mappingMode::FUNCTION1:
        {
            scalar x1Dash = mappingValuePtr_->value(x1);
            scalar x2Dash = mappingValuePtr_->value(x2);

            return value_->integrate(x1Dash, x2Dash);
        }
        case mappingMode::MINMAX:
        {
            scalar xlim0 = clamp(x1, min_, max_);
            scalar xlim1 = clamp(x2, min_, max_);

            Type intValue = value_->integrate(xlim0, xlim1);

            if (x1 < min_)
            {
                intValue += (min(min_, x2) - x1)*this->value(min_);
            }

            if (x2 > max_)
            {
                intValue += (x2 - max(max_, x1))*this->value(max_);
            }

            return intValue;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << mappingModeNames_[mappingMode_]
                << ".  Available options are: " << mappingModeNames_.sortedToc()
                << abort(FatalError);
        }
    }

    return pTraits<Type>::zero;
}


// ************************************************************************* //
