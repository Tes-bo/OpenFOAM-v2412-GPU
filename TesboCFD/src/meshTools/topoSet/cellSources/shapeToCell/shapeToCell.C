/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "shapeToCell.H"
#include "polyMesh.H"
#include "unitConversion.H"
#include "hexMatcher.H"
#include "cellFeatures.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shapeToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, shapeToCell, word);
    addToRunTimeSelectionTable(topoSetSource, shapeToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, shapeToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, shapeToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::shapeToCell::usage_
(
    shapeToCell::typeName,
    "\n    Usage: shapeToCell tet|pyr|prism|hex|tetWedge|wedge|splitHex\n\n"
    "    Select all cells of given cellShape.\n"
    "    (splitHex hardcoded with internal angle < 10 degrees)\n"
);


// Angle for polys to be considered splitHexes.
Foam::scalar Foam::shapeToCell::featureCos = Foam::cos(degToRad(10.0));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shapeToCell::combine(topoSet& set, const bool add) const
{
    if (shape_ == "splitHex")
    {
        for (label celli = 0; celli < mesh_.nCells(); ++celli)
        {
            cellFeatures superCell(mesh_, featureCos, celli);

            if (hexMatcher::test(superCell.faces()))
            {
                addOrDelete(set, celli, add);
            }
        }
    }
    else
    {
        const cellModel& wantedModel = cellModel::ref(shape_);

        const cellShapeList& cellShapes = mesh_.cellShapes();

        forAll(cellShapes, celli)
        {
            if (cellShapes[celli].model() == wantedModel)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapeToCell::shapeToCell
(
    const polyMesh& mesh,
    const word& shapeName
)
:
    topoSetCellSource(mesh),
    shape_(shapeName)
{
    if (!cellModel::ptr(shape_) && shape_ != "splitHex")
    {
        FatalErrorInFunction
            << "Illegal cell shape " << shape_ << exit(FatalError);
    }
}


Foam::shapeToCell::shapeToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    shape_(dict.getCompat<word>("shape", {{"type", 1806}}))
{}


Foam::shapeToCell::shapeToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    shape_(checkIs(is))
{
    if (!cellModel::ptr(shape_) && shape_ != "splitHex")
    {
        FatalErrorInFunction
            << "Illegal cell shape " << shape_ << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shapeToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all " << shape_ << " cells ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all " << shape_ << " cells ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
