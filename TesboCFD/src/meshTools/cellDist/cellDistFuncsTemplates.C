/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016,2024 OpenFOAM Foundation
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

#include "cellDistFuncs.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::labelHashSet Foam::cellDistFuncs::getPatchIDs() const
{
    const polyBoundaryMesh& bMesh = mesh().boundaryMesh();

    labelHashSet patchIDs(bMesh.size());

    forAll(bMesh, patchi)
    {
        if (isA<Type>(bMesh[patchi]))
        {
            patchIDs.insert(patchi);
        }
    }
    return patchIDs;
}


template<class PatchType>
Foam::scalar Foam::cellDistFuncs::smallestDist
(
    const point& p,
    const PatchType& patch,
    const labelUList& wallFaces,
    label& minFacei
) const
{
    // Return smallest true distance from p to any of wallFaces.
    // Note that even if normal hits face we still check other faces.

    const pointField& points = patch.points();

    scalar minDist = GREAT;
    minFacei = -1;

    for (const label patchFacei : wallFaces)
    {
        const pointHit curHit = patch[patchFacei].nearestPoint(p, points);

        if (curHit.distance() < minDist)
        {
            minDist = curHit.distance();
            minFacei = patchFacei;
        }
    }

    return minDist;
}


template<class PatchType>
void Foam::cellDistFuncs::getPointNeighbours
(
    const PatchType& patch,
    const label patchFacei,
    DynamicList<label>& neighbours
) const
{
    // Get point neighbours of facei (including facei). Returns number of faces.
    // Note: does not allocate storage but does use linear search to determine
    // uniqueness. For polygonal faces this might be quite inefficient.

    neighbours.clear();

    // Add myself
    neighbours.append(patchFacei);

    // Add all face neighbours
    const labelList& faceNeighbours = patch.faceFaces()[patchFacei];

    for (const label nbr : faceNeighbours)
    {
        neighbours.push_uniq(nbr);
    }

    // Add all point-only neighbours by linear searching in edge neighbours.
    // Assumes that point-only neighbours are not using multiple points on
    // face.

    const face& f = patch.localFaces()[patchFacei];

    forAll(f, fp)
    {
        label pointi = f[fp];

        const labelList& pointNbs = patch.pointFaces()[pointi];

        for (const label facei : pointNbs)
        {
            // Check for facei in edge-neighbours part of neighbours
            neighbours.push_uniq(facei);
        }
    }


    if (debug)
    {
        // Check for duplicates

        // Use hashSet to determine nbs.
        labelHashSet nbs(4*f.size());

        forAll(f, fp)
        {
            const labelList& pointNbs = patch.pointFaces()[f[fp]];
            nbs.insert(pointNbs);
        }

        // Subtract ours.
        for (const label nb : neighbours)
        {
            if (!nbs.found(nb))
            {
                SeriousErrorInFunction
                    << "getPointNeighbours : patchFacei:" << patchFacei
                    << " verts:" << f << endl;

                forAll(f, fp)
                {
                    SeriousErrorInFunction
                        << "point:" << f[fp] << " pointFaces:"
                        << patch.pointFaces()[f[fp]] << endl;
                }

                for (const label facei : neighbours)
                {
                    SeriousErrorInFunction
                        << "fast nbr:" << facei
                        << endl;
                }

                FatalErrorInFunction
                    << "Problem: fast pointNeighbours routine included " << nb
                    << " which is not in proper neighbour list " << nbs.toc()
                    << abort(FatalError);
            }
            nbs.erase(nb);
        }

        if (nbs.size())
        {
            FatalErrorInFunction
                << "Problem: fast pointNeighbours routine did not find "
                << nbs.toc() << abort(FatalError);
        }
    }
}


// ************************************************************************* //
