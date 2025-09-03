/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

Class
    Foam::polyMeshTools

Description
    Collection of static functions operating on pointMesh.

SourceFiles
    pointMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointMeshTools_H
#define Foam_pointMeshTools_H

#include "pointMesh.H"
#include "pointConstraint.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class pointMeshTools Declaration
\*---------------------------------------------------------------------------*/

class pointMeshTools
{
public:

    //- Analyse patch for feature edges, feature points. Handles points
    //- not being on a face of patch but coupled to it.
    static void featurePointsEdges
    (
        const polyMesh& mesh,

        const uindirectPrimitivePatch& boundary,
        // Per boundary face to zone
        const labelUList& faceToZone,
        // Number of zones
        const label nZones,

        const scalar edgeFeatureAngle,
        //const scalar pointFeatureAngle, //not yet done

        // Feature edge(points) internal to a zone
        labelListList& zoneToMeshPoints,
        List<pointConstraintList>& zoneToConstraints,

        // Feature edge(points) in between zones
        labelList& twoZoneMeshPoints,
        pointConstraintList& twoZoneConstraints,

        // Feature points on > 2 zones
        labelList& multiZoneMeshPoints,
        pointConstraintList& multiZoneConstraints
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
