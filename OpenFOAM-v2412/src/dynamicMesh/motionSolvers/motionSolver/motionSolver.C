/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "motionSolver.H"
#include "Time.H"
#include "polyMesh.H"
#include "dictionaryEntry.H"
#include "dlLibraryTable.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::motionSolver::stealRegistration
(
    const IOdictionary& dict
)
{
    IOobject io(dict);
    if (dict.registerObject())
    {
        // De-register if necessary
        const_cast<IOdictionary&>(dict).checkOut();
    }
    io.registerObject(true);

    return io;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver(const polyMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh)
{}


Foam::motionSolver::motionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    IOdictionary(stealRegistration(dict), dict),
    mesh_(mesh),
    coeffDict_(dict.optionalSubDict(type + "Coeffs"))
{}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::clone() const
{
    NotImplemented;
    return nullptr;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh,
    const IOdictionary& solverDict
)
{
    // The change from "solver" to "motionSolver" has not been
    // applied consistently. Accept both without complaint.
    const word solverName
    (
        solverDict.getCompat<word>("motionSolver", {{"solver", -1612}})
    );

    Info<< "Selecting motion solver: " << solverName << nl;

    mesh.time().libs().open
    (
        solverDict,
        "motionSolverLibs",
        dictionaryConstructorTablePtr_
    );

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "solver table is empty"
            << exit(FatalError);
    }

    auto* ctorPtr = dictionaryConstructorTable(solverName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            solverDict,
            "solver",
            solverName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<motionSolver>(ctorPtr(mesh, solverDict));
}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New(const polyMesh& mesh)
{
    IOdictionary solverDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    return New(mesh, solverDict);
}


Foam::motionSolver::iNew::iNew(const polyMesh& mesh)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::iNew::operator()
(
    Istream& is
) const
{
    dictionaryEntry dictEntry(dictionary::null, is);

    return motionSolver::New
    (
        mesh_,
        IOdictionary
        (
            IOobject
            (
                dictEntry.name() + ":meshSolver",
                mesh_.time().constant(),
                mesh_
            ),
            dictEntry
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();
    return curPoints();
}


void Foam::motionSolver::twoDCorrectPoints(pointField& p) const
{
    twoDPointCorrector::New(mesh_).correctPoints(p);
}


void Foam::motionSolver::updateMesh(const mapPolyMesh& mpm)
{}


bool Foam::motionSolver::writeObject
(
    IOstreamOption,
    const bool writeOnProc
) const
{
    return true;
}


bool Foam::motionSolver::read()
{
    if (regIOobject::read())
    {
        coeffDict_ = optionalSubDict(type() + "Coeffs");

        return true;
    }

    return false;
}


// ************************************************************************* //
