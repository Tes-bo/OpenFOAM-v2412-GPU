/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "adaptiveLinear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveLinear, 0);
    addToRunTimeSelectionTable(relaxationModel, adaptiveLinear, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveLinear::adaptiveLinear
(
    const dictionary& relaxationDict,
    const Time& runTime
)
:
    relaxationModel(typeName, relaxationDict, runTime),
    relaxationStart_(coeffDict().get<scalar>("relaxationStart")),
    relaxationEnd_(coeffDict().get<scalar>("relaxationEnd")),
    lastTimeValue_(runTime_.time().timeOutputValue()),
    relaxation_(relaxationStart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::adaptiveLinear::relaxation()
{
    if (runTime_.time().timeOutputValue() > lastTimeValue_)
    {
        scalar currentRelaxation = relaxation_;

        relaxation_ -=
            (relaxation_ - relaxationEnd_)
           /(
                (
                    runTime_.time().endTime().value()
                  - runTime_.time().timeOutputValue()
                )
               /(runTime_.time().timeOutputValue() - lastTimeValue_)
              + 1
            );

        lastTimeValue_ = runTime_.time().timeOutputValue();

        return currentRelaxation;
    }

    return relaxation_;
}


// ************************************************************************* //
