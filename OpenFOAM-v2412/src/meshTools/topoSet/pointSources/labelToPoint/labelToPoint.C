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

#include "labelToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(labelToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, labelToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, labelToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, labelToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, labelToPoint, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        labelToPoint,
        word,
        label
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        labelToPoint,
        istream,
        label
    );
}


Foam::topoSetSource::addToUsageTable Foam::labelToPoint::usage_
(
    labelToPoint::typeName,
    "\n    Usage: labelToPoint (i0 i1 .. in)\n\n"
    "    Select points by label\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelToPoint::labelToPoint
(
    const polyMesh& mesh,
    const labelList& labels
)
:
    topoSetPointSource(mesh),
    labels_(labels)
{}


Foam::labelToPoint::labelToPoint
(
    const polyMesh& mesh,
    labelList&& labels
)
:
    topoSetPointSource(mesh),
    labels_(std::move(labels))
{}


Foam::labelToPoint::labelToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    labels_(dict.get<labelList>("value"))
{}


Foam::labelToPoint::labelToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    labels_(checkIs(is))
{
    check(labels_, mesh.nPoints());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding points mentioned in dictionary ..." << endl;
        }

        addOrDelete(set, labels_, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing points mentioned in dictionary ..." << endl;
        }

        addOrDelete(set, labels_, false);
    }
}


// ************************************************************************* //
