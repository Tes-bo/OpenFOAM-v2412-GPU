/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "badQualityToCell.H"
#include "polyMesh.H"
#include "motionSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(badQualityToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, badQualityToCell, word);
    addToRunTimeSelectionTable(topoSetSource, badQualityToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, badQualityToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, badQualityToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        badQualityToCell,
        word,
        badQuality
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        badQualityToCell,
        istream,
        badQuality
    );
}


Foam::topoSetSource::addToUsageTable Foam::badQualityToCell::usage_
(
    badQualityToCell::typeName,
    "\n    Usage: badQualityToCell mesh-quality-dictionary\n\n"
    "    Select all cells that do not satisfy the selection criterion\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::badQualityToCell::combine(topoSet& set, const bool add) const
{
    faceSet faces(mesh_, "meshQualityFaces", mesh_.nFaces()/100+1);
    motionSmoother::checkMesh(false, mesh_, dict_, faces);
    faces.sync(mesh_);

    for (const label facei : faces)
    {
        addOrDelete(set, mesh_.faceOwner()[facei], add);
        if (mesh_.isInternalFace(facei))
        {
            addOrDelete(set, mesh_.faceNeighbour()[facei], add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::badQualityToCell::badQualityToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh),
    dict_(dict)
{}


Foam::badQualityToCell::badQualityToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    dict_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::badQualityToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding bad-quality cells" << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing bad-quality cells" << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
