/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::triSurfaceGeoMesh

Description
    The triSurface point GeoMesh (for holding point fields).

\*---------------------------------------------------------------------------*/

#ifndef Foam_triSurfacePointGeoMesh_H
#define Foam_triSurfacePointGeoMesh_H

#include "GeoMesh.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfacePointGeoMesh
:
    public GeoMesh<triSurface>
{
public:

    // Constructors

        //- Construct from triSurface reference
        explicit triSurfacePointGeoMesh(const triSurface& mesh)
        :
            GeoMesh<triSurface>(mesh)
        {}


    // Member Functions

        //- Return size. Number of points
        static label size(const Mesh& mesh)
        {
            return mesh.points().size();
        }

        //- Return size. Number of points
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
