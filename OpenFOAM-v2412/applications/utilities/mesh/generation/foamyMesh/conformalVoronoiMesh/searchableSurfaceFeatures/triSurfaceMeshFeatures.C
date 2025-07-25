/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "triSurfaceMeshFeatures.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "surfaceFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMeshFeatures, 0);
    addToRunTimeSelectionTable
    (
        searchableSurfaceFeatures,
        triSurfaceMeshFeatures,
        dict
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMeshFeatures::triSurfaceMeshFeatures
(
    const searchableSurface& surface,
    const dictionary& dict
)
:
    searchableSurfaceFeatures(surface, dict),
    includedAngle_(dict.get<scalar>("includedAngle")),
    mode_
    (
        extendedFeatureEdgeMesh::sideVolumeTypeNames_
        [
            dict.getOrDefault<word>("meshableSide", "inside")
        ]
    )
{
    Info<< indent
        << "    Included angle = " << includedAngle_ << nl
        << "    Meshable region = "
        << extendedFeatureEdgeMesh::sideVolumeTypeNames_[mode_]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extendedFeatureEdgeMesh>
Foam::triSurfaceMeshFeatures::features() const
{
    const triSurfaceMesh& surfMesh = refCast<const triSurfaceMesh>(surface());

    surfaceFeatures sFeat(surfMesh, includedAngle_);

    // TODO: Need to read on a per region basis
    boolList surfBaffleRegions
    (
        surfMesh.patches().size(),
        (mode_ == extendedFeatureEdgeMesh::BOTH ? true : false)
    );


    return autoPtr<extendedFeatureEdgeMesh>::New
    (
        sFeat,
        surface().db(),
        surface().name() + ".extendedFeatureEdgeMesh",
        surfBaffleRegions
    );
}


// ************************************************************************* //
