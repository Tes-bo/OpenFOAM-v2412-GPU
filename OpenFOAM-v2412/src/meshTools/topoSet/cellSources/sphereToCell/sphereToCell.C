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

#include "sphereToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphereToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, word);
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, sphereToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, sphereToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        sphereToCell,
        word,
        sphere
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        sphereToCell,
        istream,
        sphere
    );
}


Foam::topoSetSource::addToUsageTable Foam::sphereToCell::usage_
(
    sphereToCell::typeName,
    "\n    Usage: sphereToCell (centreX centreY centreZ) radius\n\n"
    "    Select all cells with cellCentre within bounding sphere\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphereToCell::combine(topoSet& set, const bool add) const
{
    const tmp<pointField> tctrs(this->transform(mesh_.cellCentres()));
    const pointField& ctrs = tctrs();

    const scalar orad2 = sqr(radius_);
    const scalar irad2 = innerRadius_ > 0 ? sqr(innerRadius_) : -1;

    // Treat innerRadius == 0 like unspecified innerRadius (always accept)

    forAll(ctrs, elemi)
    {
        const scalar d2 = magSqr(ctrs[elemi] - origin_);

        if ((d2 < orad2) && (d2 > irad2))
        {
            addOrDelete(set, elemi, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const point& origin,
    const scalar radius,
    const scalar innerRadius
)
:
    topoSetCellSource(mesh),
    origin_(origin),
    radius_(radius),
    innerRadius_(innerRadius)
{}


Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    origin_(dict.getCompat<vector>("origin", {{"centre", -1806}})),
    radius_(dict.getCheck<scalar>("radius", scalarMinMax::ge(0))),
    innerRadius_
    (
        dict.getCheckOrDefault<scalar>("innerRadius", 0, scalarMinMax::ge(0))
    )
{}


Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    origin_(checkIs(is)),
    radius_(readScalar(checkIs(is))),
    innerRadius_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphereToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells within sphere,"
                << " origin = " << origin_ << ", radius = " << radius_;

            if (innerRadius_ > 0)
            {
                Info<< ", innerRadius = " << innerRadius_;
            }

            Info<< endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells within sphere,"
                << " origin = " << origin_ << ", radius = " << radius_;

            if (innerRadius_ > 0)
            {
                Info<< ", innerRadius = " << innerRadius_;
            }

            Info<< endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
