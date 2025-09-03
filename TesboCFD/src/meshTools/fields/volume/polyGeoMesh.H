/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::polyGeoMesh

Description
    The polyMesh GeoMesh for holding internal fields without an fvMesh.

\*---------------------------------------------------------------------------*/

#ifndef Foam_polyGeoMesh_H
#define Foam_polyGeoMesh_H

#include "GeoMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class polyGeoMesh
:
    public GeoMesh<polyMesh>
{
public:

    // Constructors

        //- Construct from polyMesh reference
        explicit polyGeoMesh(const polyMesh& mesh)
        :
            GeoMesh<polyMesh>(mesh)
        {}


    // Member Functions

        //- Return size. Number of cells
        static label size(const Mesh& mesh)
        {
            return mesh.nCells();
        }

        //- Return size. Number of cells
        label size() const
        {
            return size(mesh_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
