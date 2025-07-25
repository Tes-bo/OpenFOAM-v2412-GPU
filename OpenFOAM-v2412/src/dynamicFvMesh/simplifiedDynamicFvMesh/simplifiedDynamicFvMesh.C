/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "simplifiedDynamicFvMesh.H"
#include "staticFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedMeshes
{
defineTypeNameAndDebug(simplifiedDynamicFvMeshBase, 0);
defineRunTimeSelectionTable(simplifiedDynamicFvMeshBase, time);
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh>
Foam::simplifiedMeshes::simplifiedDynamicFvMeshBase::New
(
    const IOobject& io
)
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        io.time().constant(),
        polyMesh::regionName(io.name()),
        io.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictHeader);

        const word modelType(dict.get<word>("dynamicFvMesh"));

        auto* ctorPtr = timeConstructorTable(modelType);

        if (ctorPtr)
        {
            Info<< "Selecting simplified mesh model " << modelType << endl;
            return autoPtr<dynamicFvMesh>(ctorPtr(io.time(), io.name()));
        }
    }

    Info<< "Selecting simplified mesh model " << staticFvMesh::typeName << endl;
    return autoPtr<dynamicFvMesh>
    (
        new SimplifiedDynamicFvMesh<staticFvMesh>(io.time(), io.name())
    );
}


// ************************************************************************* //
