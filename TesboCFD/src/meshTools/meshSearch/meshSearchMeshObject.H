/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::meshSearchMeshObject

Description
    MeshObject wrapper around meshSearch(mesh).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSearchMeshObject_H
#define meshSearchMeshObject_H

#include "MeshObject.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class meshSearchMeshObject Declaration
\*---------------------------------------------------------------------------*/

class meshSearchMeshObject
:
    public MeshObject<polyMesh, GeometricMeshObject, meshSearchMeshObject>,
    public meshSearch
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            GeometricMeshObject,
            meshSearchMeshObject
        > MeshObject_type;

public:

    // Declare name of the class and its debug switch
    TypeName("meshSearchMeshObject");


    // Constructors

        //- Construct given polyMesh
        explicit meshSearchMeshObject(const polyMesh& mesh)
        :
            MeshObject_type(mesh),
            meshSearch(mesh)
        {}


    //- Destructor
    virtual ~meshSearchMeshObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
