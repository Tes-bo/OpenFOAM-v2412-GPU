/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "cubeRootVolDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(cubeRootVolDelta, 0);
    addToRunTimeSelectionTable(LESdelta, cubeRootVolDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::cubeRootVolDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    label nD = mesh.nGeometricD();

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = deltaCoeff_ * cbrt(mesh.V());
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        const Vector<label>& directions = mesh.geometricD();

        scalar thickness = 0.0;
        for (direction dir=0; dir<directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh.bounds().span()[dir];
                break;
            }
        }

        delta_.primitiveFieldRef() = deltaCoeff_*sqrt(mesh.V()/thickness);
    }
    else
    {
        if (debug)
        {
            FatalErrorInFunction
                << "Case is not 3D or 2D, LES is not applicable"
                << exit(FatalError);
        }
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::cubeRootVolDelta::cubeRootVolDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::cubeRootVolDelta::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "deltaCoeff",
        deltaCoeff_
    );

    calcDelta();
}


void Foam::LESModels::cubeRootVolDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
