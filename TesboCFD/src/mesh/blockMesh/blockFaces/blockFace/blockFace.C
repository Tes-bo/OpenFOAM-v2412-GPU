/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "blockFace.H"
#include "blockMeshTools.H"
#include "blockVertex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockFace, 0);
    defineRunTimeSelectionTable(blockFace, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFace::blockFace(const face& vertices)
:
    vertices_(vertices)
{}


Foam::blockFace::blockFace
(
    const dictionary& dict,
    const label index,
    Istream& is
)
:
    vertices_
    (
        blockMeshTools::read<label>
        (
            is,
            dict.subOrEmptyDict("namedVertices")
        )
    )
{}


Foam::autoPtr<Foam::blockFace> Foam::blockFace::clone() const
{
    NotImplemented;
    return nullptr;
}


Foam::autoPtr<Foam::blockFace> Foam::blockFace::New
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
{
    DebugInFunction << "Constructing blockFace" << endl;

    const word faceType(is);

    auto* ctorPtr = IstreamConstructorTable(faceType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "blockFace",
            faceType,
            *IstreamConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<blockFace>(ctorPtr(dict, index, geometry, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFace::write(Ostream& os, const dictionary& d) const
{
    // Write size and start delimiter
    os << vertices_.size() << token::BEGIN_LIST;

    // Write contents
    forAll(vertices_, i)
    {
        if (i) os << token::SPACE;
        blockVertex::write(os, vertices_[i], d);
    }

    // Write end delimiter
    os << token::END_LIST;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockFace& p)
{
    os << p.vertices_ << endl;

    return os;
}


// ************************************************************************* //
