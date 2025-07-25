/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "splineEdge.H"
#include "polyLine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(splineEdge, 0);

    addToRunTimeSelectionTable
    (
        blockEdge,
        splineEdge,
        Istream
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::splineEdge::splineEdge
(
    const pointField& points,
    const edge& fromTo,
    const pointField& internalPoints
)
:
    blockEdge(points, fromTo),
    CatmullRomSpline
    (
        polyLine::concat(firstPoint(), internalPoints, lastPoint())
    )
{}


Foam::blockEdges::splineEdge::splineEdge
(
    const pointField& points,
    const label from,
    const label to,
    const pointField& internalPoints
)
:
    splineEdge(points, edge(from,to), internalPoints)
{}


Foam::blockEdges::splineEdge::splineEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces&,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    CatmullRomSpline
    (
        polyLine::concat(firstPoint(), pointField(is), lastPoint())
    )
{
    token tok(is);
    is.putBack(tok);

    // Discard unused start/end tangents
    if (tok == token::BEGIN_LIST)
    {
        vector tangent0Ignored(is);
        vector tangent1Ignored(is);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::splineEdge::position(const scalar mu) const
{
    return CatmullRomSpline::position(mu);
}


Foam::scalar Foam::blockEdges::splineEdge::length() const
{
    return CatmullRomSpline::length();
}


// ************************************************************************* //
