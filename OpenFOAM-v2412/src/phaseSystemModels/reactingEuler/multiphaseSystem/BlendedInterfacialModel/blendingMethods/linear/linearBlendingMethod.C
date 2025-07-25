/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "linearBlendingMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        linear,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::linear
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict)
{
    for (const word& phaseName : phaseNames)
    {
        minFullyContinuousAlpha_.insert
        (
            phaseName,
            dimensionedScalar
            (
                IOobject::groupName("minFullyContinuousAlpha", phaseName),
                dimless,
                dict
            )
        );

        minPartlyContinuousAlpha_.insert
        (
            phaseName,
            dimensionedScalar
            (
                IOobject::groupName("minPartlyContinuousAlpha", phaseName),
                dimless,
                dict
            )
        );

        if
        (
            minFullyContinuousAlpha_[phaseName]
          < minPartlyContinuousAlpha_[phaseName]
        )
        {
            FatalErrorInFunction
                << "The supplied fully continuous volume fraction for "
                << phaseName
                << " is less than the partly continuous value."
                << endl << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase2.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase2.name()]);

    return
        min
        (
            max
            (
                (phase2 - minPartAlpha)
               /(minFullAlpha - minPartAlpha + SMALL),
                scalar(0)
            ),
            scalar(1)
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase1.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase1.name()]);

    return
        min
        (
            max
            (
                (phase1 - minPartAlpha)
               /(minFullAlpha - minPartAlpha + SMALL),
                scalar(0)
            ),
            scalar(1)
        );
}


// ************************************************************************* //
