/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::glTF::mesh

Description
    glTF mesh

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfMesh.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_mesh_H
#define foam_gltf_mesh_H

#include "foamGltfBase.H"
#include "DynamicList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    class mesh;
}

Ostream& operator<<(Ostream& os, const glTF::mesh& mesh);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                         Class glTF::mesh Declaration
\*---------------------------------------------------------------------------*/

class mesh
:
    public base
{
    // Private Data

        //- List of fields (name, accessor ID)
        DynamicList<Tuple2<string, label>> fields_;

        //- List of colours (name, accessor ID)
        DynamicList<Tuple2<string, label>> colours_;

        //- Accessor ID
        label accessorId_;


public:

    // Constructors

        //- Default construct
        mesh();

        //- Construct from name
        explicit mesh(const word& name);


    //- Destructor
    ~mesh() = default;


    // Public Member Functions

        //- Return the accessor ID
        label& accessorId() noexcept;

        //- Add a field to the mesh
        void addField(const word& name, const label accessorId);

        //- Add a colour to the mesh
        void addColour(const label accessorId);

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
