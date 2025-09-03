/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool directionalWallPointData<Type>::update
(
    const point& pt,
    const directionalWallPointData<Type>& w2,
    const scalar tol,
    TrackingData& n
)
{
    vector d(pt - w2.origin());
    d.removeCollinear(n);

    scalar dist2 = magSqr(d);

    if (this->valid(n))
    {
        vector d(pt - this->origin());
        d.removeCollinear(n);

        scalar currentDistSqr(magSqr(d));

        scalar diff = currentDistSqr - dist2;

        if (diff < 0)
        {
            // already nearer to pt
            return false;
        }

        if
        (
            (diff < SMALL)
         || ((currentDistSqr > SMALL) && (diff/currentDistSqr < tol))
        )
        {
            // don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // current not yet set so use any value
        this->distSqr() = magSqr(pt - w2.origin());
        this->origin() = w2.origin();
        this->data() = w2.data();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
template<class Type>
inline directionalWallPointData<Type>::directionalWallPointData()
:
    wallPointData<Type>()
{}


// Construct from components
template<class Type>
inline directionalWallPointData<Type>::directionalWallPointData
(
    const point& origin,
    const Type& data,
    const scalar distSqr
)
:
    wallPointData<Type>(origin, data, distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class Type>
template<class TrackingData>
inline bool directionalWallPointData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label,
    const directionalWallPointData<Type>& neighbourWallInfo,
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
inline bool directionalWallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label,
    const directionalWallPointData<Type>& neighbourWallInfo,
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
inline bool directionalWallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const directionalWallPointData<Type>& neighbourWallInfo,
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

} // End namespace Foam

// ************************************************************************* //
