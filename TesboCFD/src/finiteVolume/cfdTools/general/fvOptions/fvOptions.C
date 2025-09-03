/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "fvOptions.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(options, 0);
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::options::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Creating finite volume options from "
            << io.instance()/io.name() << nl
            << endl;

        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        // Check if the fvOptions file is in system
        io.instance() = mesh.time().system();

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Creating finite volume options from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
        }
        else
        {
            io.readOpt(IOobject::NO_READ);
        }
    }

    return io;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::options::options
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),
    optionList(mesh, *this)
{}


Foam::fv::options& Foam::fv::options::New(const fvMesh& mesh)
{
    auto* ptr = mesh.thisDb().getObjectPtr<options>(typeName);

    if (!ptr)
    {
        DebugInFunction
            << "Constructing " << typeName
            << " for region " << mesh.name() << nl;

        ptr = new options(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


bool Foam::fv::options::read()
{
    if (IOdictionary::regIOobject::read())
    {
        optionList::read(*this);
        return true;
    }

    return false;
}


// ************************************************************************* //
