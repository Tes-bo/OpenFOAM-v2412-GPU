/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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

#include "polyMeshFilterSettings.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMeshFilterSettings::polyMeshFilterSettings(const dictionary& dict)
:
    dict_(dict),
    controlMeshQuality_
    (
        dict_.getOrDefault<Switch>("controlMeshQuality", false)
    ),
    collapseEdgesCoeffDict_(dict_.subDict("collapseEdgesCoeffs")),
    collapseFacesCoeffDict_(dict_.subOrEmptyDict("collapseFacesCoeffs")),
    meshQualityCoeffDict_(dict_.subOrEmptyDict("controlMeshQualityCoeffs")),
    minLen_(collapseEdgesCoeffDict_.get<scalar>("minimumEdgeLength")),
    maxCos_
    (
        ::cos
        (
            degToRad
            (
                collapseEdgesCoeffDict_.get<scalar>("maximumMergeAngle")
            )
        )
    ),
    edgeReductionFactor_
    (
        meshQualityCoeffDict_.getOrDefault<scalar>("edgeReductionFactor", -1)
    ),
    maxIterations_
    (
        meshQualityCoeffDict_.getOrAdd<label>("maximumIterations", 1)
    ),
    maxSmoothIters_
    (
        meshQualityCoeffDict_.getOrAdd<label>("maximumSmoothingIterations", 0)
    ),
    initialFaceLengthFactor_
    (
        collapseFacesCoeffDict_.getOrAdd<scalar>("initialFaceLengthFactor", -1)
    ),
    faceReductionFactor_
    (
        meshQualityCoeffDict_.getOrAdd<scalar>("faceReductionFactor", -1)
    ),
    maxPointErrorCount_
    (
        meshQualityCoeffDict_.getOrAdd<label>("maxPointErrorCount", 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMeshFilterSettings::writeSettings(Ostream& os) const
{
    os  << "Merging:" << nl
        << "    edges with length less than " << minLen() << " metres" << nl
        << "    edges split by a point with edges in line to within "
        << radToDeg(::acos(maxCos())) << " degrees" << nl
        << "    Minimum edge length reduction factor = "
        << edgeReductionFactor() << nl
        << endl;

    if (collapseFacesCoeffDict().empty())
    {
        os  << "Face collapsing is off" << endl;
    }
    else
    {
        os  << "Face collapsing is on" << endl;
        os  << "    Initial face length factor = "<< initialFaceLengthFactor()
            << endl;
    }

    os  << "Control mesh quality = " << controlMeshQuality().c_str() << endl;

    if (controlMeshQuality())
    {
        os  << "    Minimum edge length reduction factor = "
            << edgeReductionFactor() << nl
            << "    Minimum face area reduction factor = "
            << faceReductionFactor() << endl;

        os  << "    Maximum number of collapse iterations = " << maxIterations()
            << endl;

        os  << "    Maximum number of edge/face reduction factor smoothing "
            << "iterations = " << maxSmoothIters() << endl;

        os  << "    Maximum number of times a point can contribute to bad "
            << "faces across " << nl
            << "    collapse iterations = " << maxPointErrorCount()
            << endl;
    }

    os  << "Selectively disabling wanted collapses until resulting quality"
        << " satisfies constraints in system/meshQualityDict" << nl
        << endl;
}


// ************************************************************************* //
