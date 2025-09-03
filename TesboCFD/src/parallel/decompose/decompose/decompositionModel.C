/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "decompositionModel.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionModel, 0);
}

const Foam::word Foam::decompositionModel::canonicalName("decomposeParDict");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionModel::decompositionModel
(
    const polyMesh& mesh,
    const fileName& decompDictFile,
    const dictionary* fallback
)
:
    MeshObject_type(mesh),

    IOdictionary
    (
        IOobject::selectIO
        (
            IOobject
            (
                decompositionModel::canonicalName,
                mesh.time().system(),
                mesh.local(),
                mesh.thisDb(),
                (fallback ? IOobject::READ_IF_PRESENT : IOobject::MUST_READ),
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER,
                true    //io.globalObject()
            ),
            decompDictFile
        ),
        fallback
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

const Foam::decompositionModel& Foam::decompositionModel::New
(
    const polyMesh& mesh,
    const fileName& decompDictFile,
    const dictionary* content
)
{
    return MeshObject_type::New(mesh, decompDictFile, content);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::decompositionMethod& Foam::decompositionModel::decomposer() const
{
    if (!decomposerPtr_)
    {
        decomposerPtr_ =
            decompositionMethod::New
            (
                *this,
                this->mesh().name()  // Name of mesh region
            );
    }
    return *decomposerPtr_;
}


// ************************************************************************* //
