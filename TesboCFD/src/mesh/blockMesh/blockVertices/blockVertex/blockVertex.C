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

#include "blockVertex.H"
#include "pointVertex.H"
#include "blockMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockVertex, 0);
    defineRunTimeSelectionTable(blockVertex, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blockVertex> Foam::blockVertex::clone() const
{
    NotImplemented;
    return nullptr;
}


Foam::autoPtr<Foam::blockVertex> Foam::blockVertex::New
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
{
    DebugInFunction << "Constructing blockVertex" << endl;

    token firstToken(is);

    if (firstToken.isPunctuation(token::BEGIN_LIST))
    {
        // Putback the opening bracket
        is.putBack(firstToken);

        return autoPtr<blockVertex>
        (
            new blockVertices::pointVertex(dict, index, geometry, is)
        );
    }
    else if (firstToken.isWord())
    {
        const word faceType(firstToken.wordToken());

        auto* ctorPtr = IstreamConstructorTable(faceType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "blockVertex",
                faceType,
                *IstreamConstructorTablePtr_
            ) << abort(FatalIOError);
        }

        return autoPtr<blockVertex>(ctorPtr(dict, index, geometry, is));
    }

    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <word> or '(', found "
        << firstToken.info() << nl
        << exit(FatalIOError);

    return nullptr;
}


Foam::label Foam::blockVertex::read
(
    Istream& is,
    const dictionary& dict
)
{
    const dictionary* varDictPtr = dict.findDict("namedVertices");
    if (varDictPtr)
    {
        return blockMeshTools::read(is, *varDictPtr);
    }
    return readLabel(is);
}


void Foam::blockVertex::write
(
    Ostream& os,
    const label val,
    const dictionary& d
)
{
    const dictionary* varDictPtr = d.findDict("namedVertices");
    if (varDictPtr)
    {
        blockMeshTools::write(os, val, *varDictPtr);
    }
    else
    {
        os << val;
    }
}


// ************************************************************************* //
