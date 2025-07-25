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

#include "PurePhaseModel.H"
#include "multiphaseInterSystem.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::PurePhaseModel
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName)
{
    thermoPtr_.reset
    (
        phaseThermo::New
        (
            fluid.mesh(),
            phaseName,
            basicThermo::phasePropertyName(basicThermo::dictName, phaseName)
        )
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
void Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::solveYi
(
    PtrList<Foam::volScalarField::Internal>&,
    PtrList<Foam::volScalarField::Internal>&
)
{
    NotImplemented;
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return Y_;
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return Y_;
}


template<class BasePhaseModel, class phaseThermo>
const phaseThermo& Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::
thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo& Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}


// ************************************************************************* //
