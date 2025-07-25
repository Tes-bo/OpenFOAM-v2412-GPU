/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "linear.H"
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
        maxFullyDispersedAlpha_.insert
        (
            phaseName,
            dimensionedScalar
            (
                IOobject::groupName("maxFullyDispersedAlpha", phaseName),
                dimless,
                dict
            )
        );

        maxPartlyDispersedAlpha_.insert
        (
            phaseName,
            dimensionedScalar
            (
                IOobject::groupName("maxPartlyDispersedAlpha", phaseName),
                dimless,
                dict
            )
        );

        if
        (
            maxFullyDispersedAlpha_[phaseName]
          > maxPartlyDispersedAlpha_[phaseName]
        )
        {
            FatalErrorInFunction
                << "The supplied fully dispersed volume fraction for "
                << phaseName
                << " is greater than the partly dispersed value."
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
        maxFullAlpha(maxFullyDispersedAlpha_[phase1.name()]);
    const dimensionedScalar
        maxPartAlpha(maxPartlyDispersedAlpha_[phase1.name()]);

    return
        min
        (
            max
            (
                (phase1 - maxFullAlpha)
               /(maxPartAlpha - maxFullAlpha + SMALL),
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
        maxFullAlpha(maxFullyDispersedAlpha_[phase2.name()]);
    const dimensionedScalar
        maxPartAlpha(maxPartlyDispersedAlpha_[phase2.name()]);

    return
        min
        (
            max
            (
                (maxPartAlpha - phase2)
               /(maxPartAlpha - maxFullAlpha + SMALL),
                scalar(0)
            ),
            scalar(1)
        );
}


// ************************************************************************* //
