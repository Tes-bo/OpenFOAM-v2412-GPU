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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool Foam::wallPointData<Type>::update
(
    const point& pt,
    const wallPointData<Type>& w2,
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
            // already nearer to pt
            return false;
        }

        if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
        {
            // don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // current not yet set so use any value
        distSqr() = dist2;
        origin() = w2.origin();
        data_ = w2.data();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::wallPointData<Type>::wallPointData()
:
    wallPoint(),
    data_()
{}


template<class Type>
inline Foam::wallPointData<Type>::wallPointData
(
    const point& origin,
    const Type& data,
    const scalar distSqr
)
:
    wallPoint(origin, distSqr),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool Foam::wallPointData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const vectorField& cellCentres = mesh.primitiveMesh::cellCentres();

    return update
    (
        cellCentres[thisCelli],
        neighbourWallInfo,
        tol,
        td
    );
}


// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool Foam::wallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool Foam::wallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// ************************************************************************* //
