/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

#include "regularisationRadiusIsotropic.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isotropic, 1);
    addToRunTimeSelectionTable(regularisationRadius, isotropic, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::isotropic::computeRadius(const dictionary& dict)
{
    scalar averageVol(gAverage(mesh_.V().field()));
    const Vector<label>& geometricD = mesh_.geometricD();
    const boundBox& bounds = mesh_.bounds();
    forAll(geometricD, iDir)
    {
        if (geometricD[iDir] == -1)
        {
            averageVol /= bounds.span()[iDir];
        }
    }
    scalar radius = pow(averageVol, scalar(1)/scalar(mesh_.nGeometricD()));

    scalar multMeanRadius = dict.getOrDefault<scalar>("meanRadiusMult", 10);
    Info<< "Computed a mean radius of " << radius
        << " and multiplying with " << multMeanRadius << endl;
    return multMeanRadius*radius;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isotropic::isotropic
(
    const fvMesh& mesh,
    const dictionary& dict,
    bool adjustWallThickness
)
:
    regularisationRadius(mesh, dict, adjustWallThickness),
    radius_
    (
        "radius",
        dimLength,
        scalar
        (
            dict_.getOrDefault<scalar>("radius", computeRadius(dict))
           /(2.*::sqrt(3.))
        )
    )
{
    if (adjustWallThickness)
    {
        const scalar wallThicknessMult =
            dict.getOrDefault<scalar>("wallThicknessMult", 0.75);
        DebugInfo<<
            "Adjusting wall thickness by " << wallThicknessMult << endl;
        radius_ *= wallThicknessMult;
    }
    DebugInfo
        << "Using radius " << radius_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isotropic::addRegularisationTerm
(
    fvScalarMatrix& matrix,
    bool isTopoField
) const
{
    const volScalarField& field = matrix.psi();
    matrix -= fvm::laplacian(sqr(radius_), field);
}


// ************************************************************************* //
