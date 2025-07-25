/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::Lee
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", inv(dimTime), dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    alphaMin_(dict.getOrDefault<scalar>("alphaMin", 0))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::Kexp
(
    const volScalarField& refValue
)
{
    {
        const volScalarField from
        (
            clamp(this->pair().from(), zero_one{})
        );

        const volScalarField coeff
        (
            C_*from*this->pair().from().rho()*pos(from - alphaMin_)
           *(refValue - Tactivate_)
           /Tactivate_
        );

        if (sign(C_.value()) > 0)
        {
            return
            (
                coeff*pos(refValue - Tactivate_)
            );
        }
        else
        {
            return
            (
                coeff*pos(Tactivate_ - refValue)
            );
        }
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField from
        (
            clamp(this->pair().from(), zero_one{})
        );

        const volScalarField coeff
        (
            C_*from*this->pair().from().rho()*pos(from - alphaMin_)
            /Tactivate_
        );

        if (sign(C_.value()) > 0)
        {
            return
            (
                coeff*pos(refValue - Tactivate_)
            );
        }
        else
        {
            return
            (
                -coeff*pos(Tactivate_ - refValue)
            );
        }
    }
    else
    {
        return nullptr;
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField from
        (
            clamp(this->pair().from(), zero_one{})
        );

        const volScalarField coeff
        (
            C_*from*this->pair().from().rho()*pos(from - alphaMin_)
        );

        if (sign(C_.value()) > 0)
        {
            return
            (
                -coeff*pos(refValue - Tactivate_)
            );
        }
        else
        {
            return
            (
                -coeff*pos(Tactivate_ - refValue)
            );
        }
    }
    else
    {
        return nullptr;
    }
}


// ************************************************************************* //
