/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::simplifiedFvMesh

Description
    Functions to generate simplified finite volume meshes

SourceFiles
    simplifiedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedFvMesh_H
#define simplifiedFvMesh_H

#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simplifiedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class simplifiedFvMesh
:
    public fvMesh
{
public:

    //- Runtime type information
    TypeName("simplifiedFvMesh");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            simplifiedFvMesh,
            time,
            (
                const Time& runTime
            ),
            (runTime)
        );


    // Selectors

        //- Return a reference to the selected simplified mesh
        static autoPtr<simplifiedFvMesh> New
        (
            const word& modelType,
            const Time& runTime
        );


    //- Constructor
    simplifiedFvMesh
    (
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        labelList&& allOwner,
        labelList&& allNeighbour
    );


    // Member Functions

        //- Helper function to see if the patch type exists in the run-time
        //- selection tables
        static bool fvPatchFieldExists(const word& patchType);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
