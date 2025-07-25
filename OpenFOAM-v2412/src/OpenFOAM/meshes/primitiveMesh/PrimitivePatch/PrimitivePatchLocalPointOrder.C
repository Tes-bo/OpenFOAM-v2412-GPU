/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Description
     Orders the local points on the patch for most efficient search

\*---------------------------------------------------------------------------*/

#include "boolList.H"
#include "CircularBuffer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcLocalPointOrder() const
{
    // Note: Cannot use bandCompressing as point-point addressing does
    // not exist and is not considered generally useful.

    DebugInFunction << "Calculating local point order" << endl;

    if (localPointOrderPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "local point order already calculated"
            << abort(FatalError);
    }

    const List<face_type>& lf = localFaces();

    const labelListList& ff = faceFaces();

    localPointOrderPtr_.reset(new labelList(meshPoints().size(), -1));
    auto& pointOrder = *localPointOrderPtr_;

    boolList visitedFace(lf.size(), false);
    boolList visitedPoint(pointOrder.size(), false);

    label nPoints = 0;

    // FIFO buffer managing point/face insertion order
    CircularBuffer<label> faceOrder(32);

    forAll(lf, facei)
    {
        if (!visitedFace[facei])
        {
            faceOrder.push_back(facei);

            while (!faceOrder.empty())
            {
                // Process as FIFO
                const label curFace = faceOrder.front();
                faceOrder.pop_front();

                if (!visitedFace[curFace])
                {
                    visitedFace[curFace] = true;

                    // mark points
                    for (const label pointi : lf[curFace])
                    {
                        if (!visitedPoint[pointi])
                        {
                            visitedPoint[pointi] = true;

                            pointOrder[nPoints] = pointi;

                            ++nPoints;
                        }
                    }

                    // Add unvisited face neighbours to the list

                    for (const label nbrFacei : ff[curFace])
                    {
                        if (!visitedFace[nbrFacei])
                        {
                            faceOrder.push_back(nbrFacei);
                        }
                    }
                }
            }
        }
    }

    DebugInfo << "Calculated local point order" << endl;
}


// ************************************************************************* //
