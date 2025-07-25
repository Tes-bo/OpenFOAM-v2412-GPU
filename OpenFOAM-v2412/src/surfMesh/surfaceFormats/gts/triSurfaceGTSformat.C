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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeGTS
(
    const fileName& filename,
    const bool sort
) const
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    // Write header
    os  << "# GTS file" << endl
        << "# Regions:" << endl;

    labelList faceMap;
    surfacePatchList patches(calcPatches(faceMap));

    // Print patch names as comment
    forAll(patches, patchi)
    {
        os  << "#     " << patchi << "    "
            << patches[patchi].name() << nl;
    }
    os  << "#" << nl;

    const pointField& pts = points();

    os  << "# nPoints  nEdges  nTriangles" << nl
        << pts.size() << ' ' << nEdges() << ' ' << size() << nl;

    // Write vertex coords
    for (const point& pt : pts)
    {
        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = edges();
    const labelList& meshPts = meshPoints();

    for (const edge& e : es)
    {
        os  << meshPts[e.start()] + 1 << ' '
            << meshPts[e.end()] + 1 << nl;
    }

    // Write faces in terms of edges.
    const labelListList& faceEs = faceEdges();

    if (sort)
    {
        label faceIndex = 0;
        for (const surfacePatch& p : patches)
        {
            for (label nLocal = p.size(); nLocal--; ++faceIndex)
            {
                const label facei = faceMap[faceIndex];

                const labelList& fEdges = faceEdges()[facei];

                os  << fEdges[0] + 1 << ' '
                    << fEdges[1] + 1 << ' '
                    << fEdges[2] + 1 << ' '
                    << (*this)[facei].region() << nl;
            }
        }
    }
    else
    {
        forAll(faceEs, facei)
        {
            const labelList& fEdges = faceEdges()[facei];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << (*this)[facei].region() << nl;
        }
    }
}


// ************************************************************************* //
