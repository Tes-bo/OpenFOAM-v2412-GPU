/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "bezier.H"
#include "polyLine.H"
#include "SubList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(bezier, 0);
    addToRunTimeSelectionTable(blockEdge, bezier, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::bezier::bezier
(
    const pointField& points,
    const edge& fromTo,
    const pointField& control
)
:
    blockEdge(points, fromTo),
    control_(control)
{}


Foam::blockEdges::bezier::bezier
(
    const pointField& points,
    const label from,
    const label to,
    const pointField& control
)
:
    bezier(points, edge(from, to), control)
{}


Foam::blockEdges::bezier::bezier
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces&,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    control_
    (
        polyLine::concat(firstPoint(), pointField(is), lastPoint())
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::bezier::position(const scalar lambda) const
{
    pointField working(control_);

    label nWorking(working.size());

    forAll(working, workingI)
    {
        --nWorking;

        SubList<point>(working, nWorking) =
            (1 - lambda)*SubList<point>(working, nWorking)
          + lambda*SubList<point>(working, nWorking, 1);
    }

    return working[0];
}


Foam::scalar Foam::blockEdges::bezier::length() const
{
    NotImplemented;
    return 1;
}


// ************************************************************************* //
