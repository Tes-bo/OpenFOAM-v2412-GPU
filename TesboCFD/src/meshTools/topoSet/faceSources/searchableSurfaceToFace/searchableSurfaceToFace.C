/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "searchableSurfaceToFace.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToFace, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToFace,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetFaceSource,
        searchableSurfaceToFace,
        word
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        searchableSurfaceToFace,
        word,
        surface
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToFace::usage_
(
    searchableSurfaceToFace::typeName,
    "\n    Usage: searchableSurfaceToFace surface\n\n"
    "    Select faces with centres enclosed by the surface"
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
    word surfaceName
)
{
    dict.readIfPresent("surfaceName", surfaceName);

    if (surfaceName.empty()) surfaceName = "surface";  // failsafe

    return surfaceName;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceToFace::combine(topoSet& set, const bool add) const
{
    if (!surf_)
    {
        return;
    }
    const tmp<pointField> tctrs(this->transform(mesh_.faceCentres()));
    const pointField& ctrs = tctrs();

    const searchableSurface& s = *surf_;

    // Face centres within the enclosing volumes

    List<volumeType> volTypes;
    s.getVolumeType(ctrs, volTypes);

    const label len = volTypes.size();
    for (label elemi = 0; elemi < len; ++elemi)
    {
        if (volTypes[elemi] == volumeType::INSIDE)
        {
            addOrDelete(set, elemi, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToFace::searchableSurfaceToFace
(
    const word& surfaceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    surf_
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
{
    // Check/warn for non-enclosed
    if (surf_ && !surf_->hasVolumeType())
    {
        WarningInFunction
            << nl << "The surface " << surf_->name() << " (type: "
            << surf_->type() << ") appears to be unclosed ... ignoring"
            << nl << endl;

        surf_.clear();
    }
}


Foam::searchableSurfaceToFace::searchableSurfaceToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    searchableSurfaceToFace
    (
        dict.getCompat<word>("surfaceType", {{"surface", 0}}),
        mesh,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!surf_ || !surf_->hasVolumeType())
    {
        return;
    }

    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces enclosed by surface '"
                << surf_->name() << "' (type: " << surf_->type() << ") ..."
                << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces enclosed by surface '"
                << surf_->name() << "' (type: " << surf_->type() << ") ..."
                << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
