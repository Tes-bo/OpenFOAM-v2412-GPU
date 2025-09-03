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

Class
    Foam::petscControls

Description
    Controls for PETSc as well as global initialize/finalize.
    Uses the "system/petscControlDict" and "system/petscOptions" files
    when present.

SourceFiles
    petscControls.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscControls_H
#define Foam_module_petscControls_H

#include "IOdictionary.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                        Class petscControls Declaration
\*---------------------------------------------------------------------------*/

class petscControls
:
    public MeshObject<Time, TopologicalMeshObject, petscControls>,
    public IOdictionary
{
    // Private Data

        //- Track global state. (0=false, -1=external, +1=ours)
        static int loaded_;


public:

    //- Runtime information
    TypeName("petscControlDict");


    // Constructors

        //- Construct on top-level Time with system dictionary name.
        //  Invokes start()
        explicit petscControls(const Time& runTime);


    //- Destructor - invokes stop()
    virtual ~petscControls();


    // Member Functions

        //- Singleton access to petscControlDict on Time
        static const petscControls& New(const Time& runTime);

        //- Call PetscInitialize() if required
        static void start(const fileName& optionsFile = "");

        //- Call PetscFinalize() if required
        static void stop();

        //- Test if loaded
        bool valid() const
        {
            return loaded_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
