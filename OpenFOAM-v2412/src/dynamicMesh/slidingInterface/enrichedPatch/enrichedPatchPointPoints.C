/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "enrichedPatch.H"
#include "primitiveMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcPointPoints() const
{
    // Calculate point-point addressing
    if (pointPointsPtr_)
    {
        FatalErrorInFunction
            << "Point-point addressing already calculated."
            << abort(FatalError);
    }

    // Algorithm:
    // Go through all faces and add the previous and next point as the
    // neighbour for each point. While inserting points, reject the
    // duplicates (as every internal edge will be visited twice).
    List<DynamicList<label>> pp(meshPoints().size());

    const faceList& lf = localFaces();

    for (const face& curFace : lf)
    {
        forAll(curFace, pointi)
        {
            DynamicList<label>& curPp = pp[curFace[pointi]];

            // Do next label
            curPp.push_uniq(curFace.nextLabel(pointi));

            // Do previous label
            curPp.push_uniq(curFace.prevLabel(pointi));
        }
    }

    // Re-pack the list
    pointPointsPtr_.reset(new labelListList(pp.size()));
    auto& ppAddr = *pointPointsPtr_;

    forAll(pp, pointi)
    {
        ppAddr[pointi].transfer(pp[pointi]);
    }
}


// ************************************************************************* //
