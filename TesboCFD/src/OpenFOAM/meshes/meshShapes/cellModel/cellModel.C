/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "cellModel.H"
#include "pyramid.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::cellModel::centre
(
    const labelList& pointLabels,
    const UList<point>& points
) const
{
    // Estimate cell centre by averaging the cell points
    vector cEst = Zero;
    for (const label pointi : pointLabels)
    {
        cEst += points[pointi];
    }
    cEst /= scalar(pointLabels.size());


    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres

    scalar sumV = 0;
    vector sumVc = Zero;

    forAll(faces_, facei)
    {
        const Foam::face f(pointLabels, faces_[facei]);

        const scalar pyrVol = pyramidPointFaceRef(f, cEst).mag(points);

        if (pyrVol > SMALL)
        {
            WarningInFunction
                << "zero or negative pyramid volume: " << -pyrVol
                << " for face " << facei
                << endl;
        }

        sumV -= pyrVol;
        sumVc -= pyrVol * pyramidPointFaceRef(f, cEst).centre(points);
    }

    return sumVc/(sumV + VSMALL);
}


Foam::scalar Foam::cellModel::mag
(
    const labelList& pointLabels,
    const UList<point>& points
) const
{
    // Estimate cell centre by averaging the cell points
    vector cEst = Zero;
    for (const label pointi : pointLabels)
    {
        cEst += points[pointi];
    }
    cEst /= scalar(pointLabels.size());


    // Calculate the magnitude by summing the -mags of the pyramids
    // The sign change is because the faces point outwards
    // and a pyramid is constructed from an inward pointing face
    // and the base centre-apex vector

    scalar sumV = 0;

    forAll(faces_, facei)
    {
        const Foam::face f(pointLabels, faces_[facei]);

        const scalar pyrVol = pyramidPointFaceRef(f, cEst).mag(points);

        if (pyrVol > SMALL)
        {
            WarningInFunction
                << "zero or negative pyramid volume: " << -pyrVol
                << " for face " << facei
                << endl;
        }

        sumV -= pyrVol;
    }

    return sumV;
}


// ************************************************************************* //
