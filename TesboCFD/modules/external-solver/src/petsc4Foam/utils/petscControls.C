/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "Time.H"
#include "OSspecific.H"

#include "petscControls.H"
#include "petscErrorHandling.H"
#include "petscsys.h"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(petscControls, 0);
}

int Foam::petscControls::loaded_ = 0;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::petscControls::start(const fileName& optionsFile)
{
    int err = 0;

    if (!loaded_)
    {
        PetscBool called = PETSC_FALSE;
        err = PetscInitialized(&called);

        if (called == PETSC_TRUE)
        {
            // Someone else already called it
            loaded_ = -1;
            Info<< "PETSc already initialized - inserting OpenFOAM options"
                << nl;
            AssertPETSc(PetscOptionsInsertFile
            (
                PETSC_COMM_SELF,
                NULL,
                optionsFile.c_str(),
                PETSC_TRUE
            ));

        }
    }

    if (!loaded_)
    {
        Info<< "Initializing PETSc... " << flush;

        if (isFile(optionsFile))
        {
            err = PetscInitialize
            (
                NULL,
                NULL,
                optionsFile.c_str(),
                NULL
            );
        }
        else
        {
            err = PetscInitializeNoArguments();
        }

        if (!err)
        {
            Info<< "success" << nl;
            loaded_ = 1;
        }
        else
        {
            Info<< "fail" << nl;
        }
    }
    else
    {
        Info<< "PETSc already initialized" << nl;
    }
}


void Foam::petscControls::stop()
{
    if (loaded_ > 0)
    {
        Info<< "Finalizing PETSc" << nl;
        PetscFinalize();
        loaded_ = 0;
    }
    else if (!loaded_)
    {
        Info<< "PETSc already finalized" << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::petscControls::petscControls(const Time& runTime)
:
    MeshObject<Time, TopologicalMeshObject, petscControls>(runTime),
    IOdictionary
    (
        IOobject
        (
            petscControls::typeName,
            runTime.system(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false // no register
        )
    )
{
    runTime.printExecutionTime(Info);
    start(runTime.system()/"petscOptions");
    runTime.printExecutionTime(Info);
}


const Foam::petscControls& Foam::petscControls::New(const Time& runTime)
{
    return MeshObject<Time, TopologicalMeshObject, petscControls>::New(runTime);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::petscControls::~petscControls()
{
    petscControls::stop();
}


// ************************************************************************* //
