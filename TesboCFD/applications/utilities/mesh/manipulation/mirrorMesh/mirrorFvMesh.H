/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::mirrorFvMesh

Description

SourceFiles
    mirrorFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mirrorFvMesh_H
#define mirrorFvMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mirrorFvMesh Declaration
\*---------------------------------------------------------------------------*/

class mirrorFvMesh
:
    public fvMesh
{
    // Private data

        //- Mirrored mesh
        autoPtr<fvMesh> mirrorMeshPtr_;

        //- Cell map
        autoPtr<labelList> cellMapPtr_;

        //- Point map
        autoPtr<labelList> pointMapPtr_;


    // Private Member Functions

        //- No copy construct
        mirrorFvMesh(const mirrorFvMesh&) = delete;

        //- No copy assignment
        void operator=(const mirrorFvMesh&) = delete;


public:

    // Constructors

        //- Construct from IOobject, using system mirrorMeshDict
        explicit mirrorFvMesh(const IOobject& io);

        //- Construct from IOobject and specified mirrorMeshDict
        mirrorFvMesh(const IOobject& io, const IOdictionary& mirrorDict);


    //- Destructor
    ~mirrorFvMesh() = default;


    // Member Functions

        //- Return reference to mirror mesh
        const fvMesh& mirrorMesh() const
        {
            return mirrorMeshPtr_();
        }

        //- Mirrorred to original cell
        const labelList& cellMap() const
        {
            return cellMapPtr_();
        }

        //- Mirrorred to original point
        const labelList& pointMap() const
        {
            return pointMapPtr_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
