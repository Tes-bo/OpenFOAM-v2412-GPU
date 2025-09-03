/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::cellAspectRatio

Description
    (Rough approximation of) cell aspect ratio

SourceFiles
    cellAspectRatio.C

\*---------------------------------------------------------------------------*/

#ifndef cellAspectRatio_H
#define cellAspectRatio_H

#include "MeshObject.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cellAspectRatio Declaration
\*---------------------------------------------------------------------------*/

class cellAspectRatio
:
    public MeshObject<polyMesh, MoveableMeshObject, cellAspectRatio>,
    public scalarField
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            MoveableMeshObject,
            cellAspectRatio
        > MeshObject_type;


    // Private Member Functions

        //- Construct aspect ratio
        void calcAspectRatio();


public:

    // Declare name of the class and its debug switch
    TypeName("cellAspectRatio");


    // Constructors

        //- Construct given an polyMesh
        explicit cellAspectRatio(const polyMesh&);


    //- Destructor
    virtual ~cellAspectRatio();


    // Member functions

        //- Ignore mesh motion for now
        virtual bool movePoints()
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
