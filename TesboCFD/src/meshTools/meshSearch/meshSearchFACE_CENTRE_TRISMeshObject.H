/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::meshSearchFACE_CENTRE_TRISMeshObject

Description
    MeshObject wrapper around meshSearch(mesh,  polyMesh::FACE_CENTRE_TRIS).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSearchFACE_CENTRE_TRISMeshObject_H
#define meshSearchFACE_CENTRE_TRISMeshObject_H

#include "MeshObject.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class meshSearchFACE_CENTRE_TRISMeshObject Declaration
\*---------------------------------------------------------------------------*/

class meshSearchFACE_CENTRE_TRISMeshObject
:
    public MeshObject
    <
        polyMesh,
        GeometricMeshObject,
        meshSearchFACE_CENTRE_TRISMeshObject
    >,
    public meshSearch
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            GeometricMeshObject,
            meshSearchFACE_CENTRE_TRISMeshObject
        > MeshObject_type;

public:

    // Declare name of the class and its debug switch
    TypeName("meshSearchFACE_CENTRE_TRISMeshObject");


    // Constructors

        //- Construct given polyMesh
        explicit meshSearchFACE_CENTRE_TRISMeshObject(const polyMesh& mesh)
        :
            MeshObject_type(mesh),
            meshSearch(mesh, polyMesh::FACE_CENTRE_TRIS)
        {}


    //- Destructor
    virtual ~meshSearchFACE_CENTRE_TRISMeshObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
