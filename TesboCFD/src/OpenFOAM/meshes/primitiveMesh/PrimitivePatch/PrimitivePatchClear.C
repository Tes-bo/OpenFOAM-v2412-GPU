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

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::clearGeom()
{
    DebugInFunction << "Clearing geometric data" << nl;

    localPointsPtr_.reset(nullptr);
    faceCentresPtr_.reset(nullptr);
    faceAreasPtr_.reset(nullptr);
    magFaceAreasPtr_.reset(nullptr);
    faceNormalsPtr_.reset(nullptr);
    pointNormalsPtr_.reset(nullptr);
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::clearTopology()
{
    DebugInFunction << "Clearing patch addressing" << nl;

    // group created and destroyed together
    if (edgesPtr_ && faceFacesPtr_ && edgeFacesPtr_ && faceEdgesPtr_)
    {
        edgesPtr_.reset(nullptr);
        faceFacesPtr_.reset(nullptr);
        edgeFacesPtr_.reset(nullptr);
        faceEdgesPtr_.reset(nullptr);
    }

    boundaryPointsPtr_.reset(nullptr);
    pointEdgesPtr_.reset(nullptr);
    pointFacesPtr_.reset(nullptr);
    edgeLoopsPtr_.reset(nullptr);
    localPointOrderPtr_.reset(nullptr);
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::clearPatchMeshAddr()
{
    DebugInFunction << "Clearing patch-mesh addressing" << nl;

    meshPointsPtr_.reset(nullptr);
    meshPointMapPtr_.reset(nullptr);
    localFacesPtr_.reset(nullptr);
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::clearOut()
{
    clearGeom();
    clearTopology();
    clearPatchMeshAddr();
}


// ************************************************************************* //
