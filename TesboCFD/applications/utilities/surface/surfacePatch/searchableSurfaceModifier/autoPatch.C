/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "autoPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFeatures.H"
#include "triSurfaceMesh.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace searchableSurfaceModifiers
{
    defineTypeNameAndDebug(autoPatch, 0);
    addToRunTimeSelectionTable
    (
        searchableSurfaceModifier,
        autoPatch,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifiers::autoPatch::autoPatch
(
    const searchableSurfaces& geometry,
    const dictionary& dict
)
:
    searchableSurfaceModifier(geometry, dict),
    featureAngle_(dict.get<scalar>("featureAngle"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::searchableSurfaceModifiers::autoPatch::modify
(
    const labelList& regions,
    searchableSurface& geom
) const
{
    triSurface& surf = refCast<triSurfaceMesh>(geom);

    surfaceFeatures set(surf, 180.0-featureAngle_);

    Info<< nl
        << "Feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;


    // Get per-edge status.
    boolList borderEdge(surf.nEdges(), false);
    forAll(set.featureEdges(), i)
    {
        borderEdge[set.featureEdges()[i]] = true;
    }


    bool changed = false;

    forAll(regions, i)
    {
        label regionI = regions[i];

        labelList subRegion(surf.size(), -1);
        label nSub = 0;

        forAll(surf, faceI)
        {
            if (surf[faceI].region() == regionI && subRegion[faceI] == -1)
            {
                PatchTools::markZone(surf, borderEdge, faceI, nSub, subRegion);
                nSub++;
            }
        }

        Info<< "For region " << surf.patches()[regionI].name()
            << " detected sub regions " << nSub << endl;


        if (nSub > 1)
        {
            label nOldPatches = surf.patches().size();

            for (label subI = 0; subI < nSub; subI++)
            {
                geometricSurfacePatch patch
                (
                    surf.patches()[regionI].name() + Foam::name(subI),
                    surf.patches().size(),
                    surf.patches()[regionI].geometricType()
                );


                Info<< "For region " << surf.patches()[regionI].name()
                    << " creating region " << patch.name() << endl;

                surf.patches().append(patch);
            }

            forAll(subRegion, faceI)
            {
                if (subRegion[faceI] != -1)
                {
                    surf[faceI].region() = nOldPatches + subRegion[faceI];
                }
            }
            changed = true;
        }
    }

    return changed;
}


// ************************************************************************* //
