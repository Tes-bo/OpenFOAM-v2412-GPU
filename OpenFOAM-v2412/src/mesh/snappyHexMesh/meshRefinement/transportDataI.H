/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::transportData::update
(
    const point& pt,
    const transportData& w2,
    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 = magSqr(pt - w2.origin());

    if (valid(td))
    {
        const scalar diff = distSqr() - dist2;

        if (diff < 0)
        {
            // Already nearer to pt
            return false;
        }

        if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
        {
            // Don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // current not yet set so use any value
        distSqr() = dist2;
        origin() = w2.origin();
        data() = w2.data();

        if (distSqr() > sqr(0.25*data()))
        {
            // No need to transport gap data since too far away
            return false;
        }

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::transportData::transportData()
:
    wallPointData<scalar>()
{}


inline Foam::transportData::transportData
(
    const point& origin,
    const scalar gapSize,
    const scalar distSqr
)
:
    wallPointData<scalar>(origin, gapSize, distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::transportData::updateCell
(
    const polyMesh& mesh,
    const label cellI,
    const label faceI,
    const transportData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (td.surfaceIndex_[faceI] != -1)
    {
        return false;
    }

    const vectorField& cellCentres = mesh.primitiveMesh::cellCentres();

    bool updated = update
    (
        cellCentres[cellI],
        neighbourWallInfo,
        tol,
        td
    );

    return updated;
}


template<class TrackingData>
inline bool Foam::transportData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const transportData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (td.surfaceIndex_[thisFaceI] != -1)
    {
        return false;
    }

    return update
    (
        mesh.faceCentres()[thisFaceI],
        neighbourWallInfo,
        tol,
        td
    );
}


template<class TrackingData>
inline bool Foam::transportData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const transportData& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (td.surfaceIndex_[thisFaceI] != -1)
    {
        return false;
    }

    return update
    (
        mesh.faceCentres()[thisFaceI],
        neighbourWallInfo,
        tol,
        td
    );
}


// ************************************************************************* //
