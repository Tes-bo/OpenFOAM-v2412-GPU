/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "displacementMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMethod, 0);
    defineRunTimeSelectionTable(displacementMethod, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMethod::displacementMethod
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    mesh_(mesh),
    patchIDs_(patchIDs),
    motionPtr_(motionSolver::New(mesh_)),
    maxDisplacement_(SMALL)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::displacementMethod> Foam::displacementMethod::New
(
    fvMesh& mesh,
    const labelList& patchIDs
)
{
    // To ensure that displacement method has the same
    // type as the motion solver, construct it based on its name
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    word solverType(dynamicMeshDict.get<word>("solver"));

    Info<< "displacementMethod type : " << solverType << endl;

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dynamicMeshDict,
            "solver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<displacementMethod>(ctorPtr(mesh, patchIDs));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::displacementMethod::preferPointField() const
{
    return true;
}


void Foam::displacementMethod::boundControlField(vectorField& controlField)
{
    // Does nothing in base
}


Foam::autoPtr<Foam::motionSolver>& Foam::displacementMethod::getMotionSolver()
{
    return motionPtr_;
}


Foam::scalar Foam::displacementMethod::getMaxDisplacement() const
{
    return maxDisplacement_;
}


void Foam::displacementMethod::update()
{
    scalar timeBef = mesh_.time().elapsedCpuTime();

    // Compute max mesh movement
    tmp<pointField> tnewPoints = motionPtr_->newPoints();
    Info<< "Max mesh movement magnitude "
        << gMax(mag(tnewPoints() - mesh_.points())) << endl;

    // Update mesh as in dynamicFvMesh
    mesh_.movePoints(tnewPoints());
    scalar timeAft = mesh_.time().elapsedCpuTime();
    Info<< "Mesh movement took " << timeAft - timeBef << " seconds" << endl;
    if (!mesh_.steady())
    {
        mesh_.moving(false);
    }
}


// ************************************************************************* //
