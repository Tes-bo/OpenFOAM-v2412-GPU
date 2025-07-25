/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "triSurface.H"
#include "mergePoints.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::stitchTriangles
(
    const scalar tol,
    bool verbose
)
{
    pointField& ps = this->storedPoints();

    // Merge points (inplace)
    labelList pointMap;
    label nChanged = Foam::inplaceMergePoints(ps, tol, verbose, pointMap);

    if (nChanged)
    {
        if (verbose)
        {
            Pout<< "stitchTriangles : Merged from " << pointMap.size()
                << " points down to " << ps.size() << endl;
        }

        // Reset the triangle point labels to the unique points array
        label newTriangleI = 0;
        forAll(*this, i)
        {
            const labelledTri& tri = operator[](i);
            labelledTri newTri
            (
                pointMap[tri[0]],
                pointMap[tri[1]],
                pointMap[tri[2]],
                tri.region()
            );

            if (newTri.good())
            {
                operator[](newTriangleI++) = newTri;
            }
            else if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removing triangle " << i
                    << " with non-unique vertices." << nl
                    << "    vertices   :" << newTri << nl
                    << "    coordinates:" << newTri.points(ps)
                    << endl;
            }
        }

        if (newTriangleI != size())
        {
            if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removed " << size() - newTriangleI
                    << " triangles" << endl;
            }
            setSize(newTriangleI);

            // And possibly compact out any unused points (since used only
            // by triangles that have just been deleted)
            // Done in two passes to save memory (pointField)

            // 1. Detect only
            bitSet pointIsUsed(ps.size());

            label nPoints = 0;

            for (const labelledTri& f : *this)
            {
                for (const label pointi : f)
                {
                    if (pointIsUsed.set(pointi))
                    {
                        ++nPoints;
                    }
                }
            }

            if (nPoints != ps.size())
            {
                // 2. Compact.
                pointMap.setSize(ps.size());
                label newPointi = 0;
                forAll(pointIsUsed, pointi)
                {
                    if (pointIsUsed[pointi])
                    {
                        ps[newPointi] = ps[pointi];
                        pointMap[pointi] = newPointi++;
                    }
                }
                ps.setSize(newPointi);

                newTriangleI = 0;
                for (const labelledTri& f : *this)
                {
                    operator[](newTriangleI++) = labelledTri
                    (
                        pointMap[f[0]],
                        pointMap[f[1]],
                        pointMap[f[2]],
                        f.region()
                    );
                }
            }
        }
    }

    return bool(nChanged);
}


// ************************************************************************* //
