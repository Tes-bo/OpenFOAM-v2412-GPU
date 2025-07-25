/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "boundaryToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, boundaryToFace, word);
    addToRunTimeSelectionTable(topoSetSource, boundaryToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, boundaryToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, boundaryToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        boundaryToFace,
        word,
        boundary
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        boundaryToFace,
        istream,
        boundary
    );
}


Foam::topoSetSource::addToUsageTable Foam::boundaryToFace::usage_
(
    boundaryToFace::typeName,
    "\n    Usage: boundaryToFace\n\n"
    "    Select all boundary faces\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryToFace::combine(topoSet& set, const bool add) const
{
    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        ++facei
    )
    {
        addOrDelete(set, facei, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryToFace::boundaryToFace(const polyMesh& mesh)
:
    topoSetFaceSource(mesh)
{}


Foam::boundaryToFace::boundaryToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict)
{}


Foam::boundaryToFace::boundaryToFace
(
    const polyMesh& mesh,
    Istream&
)
:
    topoSetFaceSource(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all boundary faces ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all boundary faces ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
