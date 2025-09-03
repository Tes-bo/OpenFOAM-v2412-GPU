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

#include "labelToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(labelToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, labelToFace, word);
    addToRunTimeSelectionTable(topoSetSource, labelToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, labelToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, labelToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        labelToFace,
        word,
        label
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        labelToFace,
        istream,
        label
    );
}


Foam::topoSetSource::addToUsageTable Foam::labelToFace::usage_
(
    labelToFace::typeName,
    "\n    Usage: labelToFace (i0 i1 .. in)\n\n"
    "    Select faces by label\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelToFace::labelToFace
(
    const polyMesh& mesh,
    const labelList& labels
)
:
    topoSetFaceSource(mesh),
    labels_(labels)
{}


Foam::labelToFace::labelToFace
(
    const polyMesh& mesh,
    labelList&& labels
)
:
    topoSetFaceSource(mesh),
    labels_(std::move(labels))
{}


Foam::labelToFace::labelToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    labels_(dict.get<labelList>("value"))
{}


Foam::labelToFace::labelToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    labels_(checkIs(is))
{
    check(labels_, mesh.nFaces());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces mentioned in dictionary"
                << " ..." << endl;
        }

        addOrDelete(set, labels_, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces mentioned dictionary"
                << " ..." << endl;
        }

        addOrDelete(set, labels_, false);
    }
}


// ************************************************************************* //
