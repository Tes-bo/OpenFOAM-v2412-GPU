/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "cvControls.H"
#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cvControls::cvControls
(
    const dictionary& foamyHexMeshDict
)
:
    foamyHexMeshDict_(foamyHexMeshDict)
{
    // Surface conformation controls

    const dictionary& surfDict
    (
        foamyHexMeshDict_.subDict("surfaceConformation")
    );

    pointPairDistanceCoeff_ =
        surfDict.get<scalar>("pointPairDistanceCoeff");

    mixedFeaturePointPPDistanceCoeff_ =
        surfDict.get<scalar>("mixedFeaturePointPPDistanceCoeff");

    featurePointExclusionDistanceCoeff_ =
        surfDict.get<scalar>("featurePointExclusionDistanceCoeff");

    featureEdgeExclusionDistanceCoeff_ =
        surfDict.get<scalar>("featureEdgeExclusionDistanceCoeff");

    surfaceSearchDistanceCoeff_ =
        surfDict.get<scalar>("surfaceSearchDistanceCoeff");

    maxSurfaceProtrusionCoeff_ =
        surfDict.get<scalar>("maxSurfaceProtrusionCoeff");

    maxQuadAngle_ = surfDict.get<scalar>("maxQuadAngle");

    surfaceConformationRebuildFrequency_ = max
    (
        1,
        surfDict.get<label>("surfaceConformationRebuildFrequency")
    );


    const dictionary& featurePointControlsDict
    (
        surfDict.subDict("featurePointControls")
    );

    specialiseFeaturePoints_ =
        featurePointControlsDict.get<Switch>("specialiseFeaturePoints");

    guardFeaturePoints_ =
        featurePointControlsDict.get<Switch>("guardFeaturePoints");

    edgeAiming_ =
        featurePointControlsDict.get<Switch>("edgeAiming");

    if (!guardFeaturePoints_)
    {
        snapFeaturePoints_ =
            featurePointControlsDict.get<Switch>("snapFeaturePoints");
    }

    circulateEdges_ =
        featurePointControlsDict.get<Switch>("circulateEdges");

    // Controls for coarse surface conformation

    const dictionary& conformationControlsDict
    (
        surfDict.subDict("conformationControls")
    );

    surfacePtExclusionDistanceCoeff_ =
        conformationControlsDict.get<scalar>("surfacePtExclusionDistanceCoeff");

    edgeSearchDistCoeffSqr_ = sqr
    (
        conformationControlsDict.get<scalar>("edgeSearchDistCoeff")
    );

    surfacePtReplaceDistCoeffSqr_ = sqr
    (
        conformationControlsDict.get<scalar>("surfacePtReplaceDistCoeff")
    );

    maxConformationIterations_ =
        conformationControlsDict.get<label>("maxIterations");

    iterationToInitialHitRatioLimit_ =
        conformationControlsDict.get<scalar>("iterationToInitialHitRatioLimit");


    // Motion control controls

    const dictionary& motionDict(foamyHexMeshDict_.subDict("motionControl"));

    defaultCellSize_ = motionDict.get<scalar>("defaultCellSize");

    minimumCellSize_ =
        motionDict.get<scalar>("minimumCellSizeCoeff")*defaultCellSize_;

    objOutput_ =
        motionDict.getOrDefault<Switch>("objOutput", false);

    timeChecks_ =
        motionDict.getOrDefault<Switch>("timeChecks", false);

    printVertexInfo_ =
        motionDict.getOrDefault<Switch>("printVertexInfo", false);

    if (Pstream::parRun())
    {
        maxLoadUnbalance_ = motionDict.get<scalar>("maxLoadUnbalance");
    }
    else
    {
        maxLoadUnbalance_ = -1;
    }

    cosAlignmentAcceptanceAngle_ = cos
    (
        degToRad(motionDict.get<scalar>("alignmentAcceptanceAngle"))
    );


    // Point removal criteria

    const dictionary& insertionDict
    (
        motionDict.subDict("pointInsertionCriteria")
    );

    insertionDistCoeff_ =
        insertionDict.get<scalar>("cellCentreDistCoeff");

    faceAreaRatioCoeff_ =
        insertionDict.get<scalar>("faceAreaRatioCoeff");

    cosInsertionAcceptanceAngle_ = cos
    (
        degToRad(insertionDict.get<scalar>("acceptanceAngle"))
    );

    // Point removal criteria

    const dictionary& removalDict
    (
        motionDict.subDict("pointRemovalCriteria")
    );

    removalDistCoeff_ =
        removalDict.get<scalar>("cellCentreDistCoeff");

    // polyMesh filtering controls

    const dictionary& filteringDict
    (
        foamyHexMeshDict_.subDict("polyMeshFiltering")
    );

    filterEdges_ =
        filteringDict.getOrDefault<Switch>("filterEdges", true);

    filterFaces_ =
        filteringDict.getOrDefault<Switch>("filterFaces", false);

    if (filterFaces_)
    {
        filterEdges_ = filterFaces_;
    }

    writeTetDualMesh_ =
        filteringDict.get<Switch>("writeTetDualMesh");

    writeCellShapeControlMesh_ =
        filteringDict.get<Switch>("writeCellShapeControlMesh");

    if (Pstream::parRun())
    {
        writeBackgroundMeshDecomposition_ =
            filteringDict.get<Switch>("writeBackgroundMeshDecomposition");
    }
    else
    {
        writeBackgroundMeshDecomposition_ = Switch::FALSE;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvControls::~cvControls()
{}


// ************************************************************************* //
