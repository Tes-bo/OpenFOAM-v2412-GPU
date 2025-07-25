/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "snappySnapDriver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList>
Foam::labelList Foam::snappySnapDriver::getFacePoints
(
    const indirectPrimitivePatch& pp,
    const FaceList& faces
)
{
    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(faces, i)
    {
        const face& f = faces[i];

        forAll(f, fp)
        {
            label meshPointi = f[fp];

            const auto iter = pp.meshPointMap().cfind(meshPointi);

            if (iter.good())
            {
                const label pointi = iter.val();
                pointOnZone[pointi] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


// ************************************************************************* //
