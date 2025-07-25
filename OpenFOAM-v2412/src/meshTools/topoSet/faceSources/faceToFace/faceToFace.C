/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "faceToFace.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, faceToFace, word);
    addToRunTimeSelectionTable(topoSetSource, faceToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, faceToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, faceToFace, istream);
}


Foam::topoSetSource::addToUsageTable Foam::faceToFace::usage_
(
    faceToFace::typeName,
    "\n    Usage: faceToFace <faceSet>\n\n"
    "    Select all faces in the faceSet\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceToFace::faceToFace
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetFaceSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false)
{}


Foam::faceToFace::faceToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_))
{}


Foam::faceToFace::faceToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all elements of face "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.addSet(mesh_.faceZones()[setName]);
            }
            else
            {
                faceSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);

                set.addSet(loadedSet);
            }
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all elements of face "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.subtractSet(mesh_.faceZones()[setName]);
            }
            else
            {
                faceSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);

                set.subtractSet(loadedSet);
            }
        }
    }
}


// ************************************************************************* //
