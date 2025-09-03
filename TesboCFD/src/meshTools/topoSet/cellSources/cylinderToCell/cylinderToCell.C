/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "cylinderToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderToCell, word);
    addToRunTimeSelectionTable(topoSetSource, cylinderToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, cylinderToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, cylinderToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        cylinderToCell,
        word,
        cylinder
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        cylinderToCell,
        istream,
        cylinder
    );
}


Foam::topoSetSource::addToUsageTable Foam::cylinderToCell::usage_
(
    cylinderToCell::typeName,
    "\n    Usage: cylinderToCell (p1X p1Y p1Z) (p2X p2Y p2Z) radius\n\n"
    "    Select all cells with cell centre within bounding cylinder\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderToCell::combine(topoSet& set, const bool add) const
{
    const tmp<pointField> tctrs(this->transform(mesh_.cellCentres()));
    const pointField& ctrs = tctrs();

    const vector axis = (point2_ - point1_);
    const scalar magAxis2 = magSqr(axis);
    const scalar orad2 = sqr(radius_);
    const scalar irad2 = innerRadius_ > 0 ? sqr(innerRadius_) : -1;

    // Treat innerRadius == 0 like unspecified innerRadius (always accept)

    forAll(ctrs, elemi)
    {
        const vector d = ctrs[elemi] - point1_;
        const scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            const scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if ((d2 < orad2) && (d2 > irad2))
            {
                addOrDelete(set, elemi, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderToCell::cylinderToCell
(
    const polyMesh& mesh,
    const point& point1,
    const point& point2,
    const scalar radius,
    const scalar innerRadius
)
:
    topoSetCellSource(mesh),
    point1_(point1),
    point2_(point2),
    radius_(radius),
    innerRadius_(innerRadius)
{
    if (innerRadius_ > radius_)
    {
        FatalErrorInFunction
            << "inner radius = " << innerRadius_ << " cannot be larger than "
            << "outer radius = " << radius_
            << exit(FatalError);
    }
}


Foam::cylinderToCell::cylinderToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    point1_(dict.getCompat<point>("point1", {{"p1", -2112}})),
    point2_(dict.getCompat<point>("point2", {{"p2", -2112}})),
    radius_(dict.getCompat<scalar>("radius", {{"outerRadius", -2112}})),
    innerRadius_
    (
        dict.getCheckOrDefault<scalar>("innerRadius", 0, scalarMinMax::ge(0))
    )
{}


Foam::cylinderToCell::cylinderToCell
(
    const polyMesh& mesh,
    Istream& is,
    const bool mandatoryInnerRadius
)
:
    topoSetCellSource(mesh),
    point1_(checkIs(is)),
    point2_(checkIs(is)),
    radius_(readScalar(checkIs(is))),
    innerRadius_(0)
{
    if (mandatoryInnerRadius)
    {
        innerRadius_ = readScalar(checkIs(is));
    }
}


Foam::cylinderToCell::cylinderToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    cylinderToCell(mesh, is, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells with centres within cylinder"
                << ", with point1 = " << point1_
                << ", point2 = " << point2_
                << ", radius = " << radius_;

            if (innerRadius_ > 0)
            {
                Info<< ", inner radius = " << innerRadius_;
            }

            Info<< endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells with centres within cylinder"
                << ", with point1 = " << point1_
                << ", point2 = " << point2_
                << ", radius = " << radius_;

            if (innerRadius_ > 0)
            {
                Info<< ", inner radius = " << innerRadius_;
            }

            Info<< endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
