/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::boundaryFaces() const
{
    labelList bndFaces(nBoundaryEdges());

    // The boundary slice
    const SubList<labelList> bndEdgeToFace
    (
        edgeFaces(),
        bndFaces.size(),
        nInternalEdges()
    );

    // By definition boundary edges have a _single_ face attached
    forAll(bndFaces, i)
    {
        bndFaces[i] = bndEdgeToFace[i][0];
    }

    return bndFaces;
}


template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::uniqBoundaryFaces() const
{
    labelList bndFaces(this->boundaryFaces());

    const label len = bndFaces.size();

    if (len > 1)
    {
        Foam::sort(bndFaces);

        label prev = bndFaces[0];
        label nUniq = 1;

        for (label i=1; i < len; ++i)
        {
            if (prev != bndFaces[i])
            {
                prev = bndFaces[i];
                bndFaces[nUniq] = prev;
                ++nUniq;
            }
        }

        bndFaces.resize(nUniq);
    }

    return bndFaces;
}


// ************************************************************************* //
