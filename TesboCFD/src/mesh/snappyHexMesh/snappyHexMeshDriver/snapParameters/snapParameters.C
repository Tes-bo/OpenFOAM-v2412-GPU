/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "snapParameters.H"
#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::snapParameters::snapParameters(const dictionary& dict, const bool dryRun)
:
    nSmoothPatch_
    (
        meshRefinement::get<label>(dict, "nSmoothPatch", dryRun)
    ),
    nSmoothInternal_(dict.getOrDefault<label>("nSmoothInternal", 0)),
    snapTol_
    (
        meshRefinement::get<scalar>(dict, "tolerance", dryRun)
    ),
    nSmoothDispl_
    (
        meshRefinement::get<label>(dict, "nSolveIter", dryRun)
    ),
    nSnap_
    (
        meshRefinement::get<label>(dict, "nRelaxIter", dryRun)
    ),
    nFeatureSnap_(dict.getOrDefault("nFeatureSnapIter", -1)),
    explicitFeatureSnap_(dict.getOrDefault("explicitFeatureSnap", true)),
    implicitFeatureSnap_(dict.getOrDefault("implicitFeatureSnap", false)),
    multiRegionFeatureSnap_
    (
        dict.getOrDefault("multiRegionFeatureSnap", false)
    ),
    detectNearSurfacesSnap_
    (
        dict.getOrDefault("detectNearSurfacesSnap", true)
    ),
    strictRegionSnap_
    (
        dict.getOrDefault("strictRegionSnap", false)
    ),
    detectBaffles_(dict.getOrDefault("detectBaffles", true)),
    baffleFeaturePoints_(dict.getOrDefault("baffleFeaturePoints", false)),
    releasePoints_(dict.getOrDefault("releasePoints", false)),
    stringFeatures_(dict.getOrDefault("stringFeatures", true)),
    avoidDiagonal_(dict.getOrDefault("avoidDiagonal", false)),
    nFaceSplitInterval_
    (
        dict.getOrDefault("nFaceSplitInterval", labelMin)
    ),
    concaveAngle_(dict.getOrDefault<scalar>("concaveAngle", 45)),
    minAreaRatio_(dict.getOrDefault<scalar>("minAreaRatio", 0.3)),
    dryRun_(dryRun)
{}


// ************************************************************************* //
