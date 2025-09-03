/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "searchableSurfaceToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "searchableSurface.H"
#include "syncTools.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToFaceZone, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToFaceZone,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        searchableSurfaceToFaceZone,
        word
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        searchableSurfaceToFaceZone,
        word,
        surface
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToFaceZone::usage_
(
    searchableSurfaceToFaceZone::typeName,
    "\n    Usage: searchableSurfaceToFaceZone surface\n\n"
    "    Select all faces whose cell-cell centre vector intersects the surface "
    "\n"
);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Difficult to get a good default name from the dictionary name.
// It could be
//     sourceInfo { .. }
// But even with something like
//     mySurf.stl { .. }
// The dictName() method will only return the "stl" ending.

static inline word getSurfaceName
(
    const dictionary& dict,
    const word& defaultName
)
{
    return
        dict.getOrDefaultCompat<word>
        (
            "surfaceName",
            {{"name", 1806}},
            defaultName
        );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToFaceZone::searchableSurfaceToFaceZone
(
    const word& surfaceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    surfacePtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                getSurfaceName(dict, mesh.objectRegistry::db().name()),
                mesh.time().constant(),     // Instance
                "triSurface",               // Local
                mesh.objectRegistry::db(),  // Registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{}


Foam::searchableSurfaceToFaceZone::searchableSurfaceToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    searchableSurfaceToFaceZone
    (
        dict.getCompat<word>("surfaceType", {{"surface", 0}}),
        mesh,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
        return;
    }
    else
    {
        faceZoneSet& fzSet = refCast<faceZoneSet>(set);

        // Get cell-cell centre vectors
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointField start(mesh_.nFaces());
        pointField end(mesh_.nFaces());

        const pointField& cc = mesh_.cellCentres();

        // Internal faces
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            start[facei] = cc[mesh_.faceOwner()[facei]];
            end[facei] = cc[mesh_.faceNeighbour()[facei]];
        }

        // Boundary faces
        vectorField nbrCellCentres;
        syncTools::swapBoundaryCellPositions(mesh_, cc, nbrCellCentres);

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        for (const polyPatch& pp : pbm)
        {
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    start[facei] = cc[mesh_.faceOwner()[facei]];
                    end[facei] = nbrCellCentres[facei-mesh_.nInternalFaces()];
                }
            }
            else
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    start[facei] = cc[mesh_.faceOwner()[facei]];
                    end[facei] = mesh_.faceCentres()[facei];
                }
            }
        }


        // Do all intersection tests
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        List<pointIndexHit> hits;
        surfacePtr_().findLine(start, end, hits);
        pointField normals;
        surfacePtr_().getNormal(hits, normals);


        // Select intersected faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        if (action == topoSetSource::ADD || action == topoSetSource::NEW)
        {
            if (verbose_)
            {
                Info<< "    Adding all faces from surface "
                    << surfacePtr_().name() << " ..." << endl;
            }

            DynamicList<label> newAddressing(fzSet.addressing());
            DynamicList<bool> newFlipMap(fzSet.flipMap());

            forAll(hits, facei)
            {
                if (hits[facei].hit() && !fzSet.found(facei))
                {
                    newAddressing.append(facei);
                    vector d = end[facei]-start[facei];
                    newFlipMap.append((normals[facei] & d) < 0);
                }
            }

            fzSet.addressing().transfer(newAddressing);
            fzSet.flipMap().transfer(newFlipMap);
            fzSet.updateSet();
        }
        else if (action == topoSetSource::SUBTRACT)
        {
            if (verbose_)
            {
                Info<< "    Removing all faces from surface "
                    << surfacePtr_().name() << " ..." << endl;
            }

            // Start off empty
            DynamicList<label> newAddressing(fzSet.addressing().size());
            DynamicList<bool> newFlipMap(fzSet.flipMap().size());

            forAll(fzSet.addressing(), i)
            {
                if (!hits[fzSet.addressing()[i]].hit())
                {
                    newAddressing.append(fzSet.addressing()[i]);
                    newFlipMap.append(fzSet.flipMap()[i]);
                }
            }
            fzSet.addressing().transfer(newAddressing);
            fzSet.flipMap().transfer(newFlipMap);
            fzSet.updateSet();
        }
    }
}


// ************************************************************************* //
