/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::staticFvMesh

Description
    Foam::staticFvMesh

SourceFiles
    staticFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticFvMesh_H
#define staticFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticFvMesh Declaration
\*---------------------------------------------------------------------------*/

class staticFvMesh
:
    public dynamicFvMesh
{
    // Private Member Functions

        //- No copy construct
        staticFvMesh(const staticFvMesh&) = delete;

        //- No copy assignment
        void operator=(const staticFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("staticFvMesh");


    // Constructors

        //- Construct from IOobject
        staticFvMesh(const IOobject& io, const bool doInit=true);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        staticFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );


    //- Destructor
    ~staticFvMesh() = default;


    // Member Functions

        //- Is mesh dynamic
        virtual bool dynamic() const
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
