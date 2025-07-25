/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
------------------------------------------------------------------------------
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

#include "log.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(log, 0);
    addToRunTimeSelectionTable(functionObject, log, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::log::calc()
{
    if (foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& x = lookupObject<volScalarField>(fieldName_);

        // Switch-off dimension checking if requested
        const bool oldDimChecking = dimensionSet::checking();

        if (!checkDimensions_)
        {
            dimensionSet::checking(false);
        }

        bool stored = store
        (
            resultName_,
            scale_*Foam::log(max(x, clipValue_)) + offset_
        );

        // Reinstate dimension checking
        if (!checkDimensions_)
        {
            dimensionSet::checking(oldDimChecking);
        }

        return stored;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::log::log
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName),
    checkDimensions_(true),
    clipValue_(SMALL),
    scale_(1.0),
    offset_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::log::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && fieldExpression::read(dict))
    {
        checkDimensions_ = dict.getOrDefault<Switch>("checkDimensions", true);
        clipValue_ =
            dict.getCheckOrDefault<scalar>
            (
                "clip",
                SMALL,
                scalarMinMax::ge(SMALL)
            );
        scale_ = dict.getOrDefault<scalar>("scale", 1.0);
        offset_ = dict.getOrDefault<scalar>("offset", 0.0);

        return true;
    }

    return false;
}


// ************************************************************************* //
