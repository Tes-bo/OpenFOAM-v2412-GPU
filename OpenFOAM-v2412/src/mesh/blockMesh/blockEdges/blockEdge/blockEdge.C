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

#include "blockEdge.H"
#include "blockVertex.H"
#include "polyLine.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockEdge, 0);
    defineRunTimeSelectionTable(blockEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdge::blockEdge
(
    const pointField& points,
    const edge& fromTo
)
:
    points_(points),
    start_(fromTo.first()),
    end_(fromTo.last())
{}


Foam::blockEdge::blockEdge
(
    const dictionary& dict,
    const label index,
    const pointField& points,
    Istream& is
)
:
    points_(points),
    start_(blockVertex::read(is, dict)),
    end_(blockVertex::read(is, dict))
{}


Foam::autoPtr<Foam::blockEdge> Foam::blockEdge::clone() const
{
    NotImplemented;
    return nullptr;
}


Foam::autoPtr<Foam::blockEdge> Foam::blockEdge::New
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
{
    DebugInFunction << "Constructing blockEdge" << endl;

    const word edgeType(is);

    auto* ctorPtr = IstreamConstructorTable(edgeType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "blockEdge",
            edgeType,
            *IstreamConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<blockEdge>(ctorPtr(dict, index, geometry, points, is));
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::pointField Foam::blockEdge::appendEndPoints
(
    const pointField& p,
    const label from,
    const label to,
    const pointField& intermediate
)
{
    return pointField(polyLine::concat(p[from], intermediate, p[to]));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::blockEdge::position(const scalarList& lambdas) const
{
    auto tpoints = tmp<pointField>::New(lambdas.size());
    auto& points = tpoints.ref();

    forAll(lambdas, i)
    {
        points[i] = position(lambdas[i]);
    }
    return tpoints;
}


void Foam::blockEdge::write(Ostream& os, const dictionary& dict) const
{
    blockVertex::write(os, start_, dict);
    os << tab;
    blockVertex::write(os, end_, dict);
    os << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockEdge& e)
{
    os << e.start_ << tab << e.end_ << endl;

    return os;
}


// ************************************************************************* //
