/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::ijkMesh

Description
    A simple i-j-k (row-major order) to linear addressing for a
    rectilinear mesh. Since the underlying mesh is rectilinear, some
    mesh-related sizing information can be derived directly from the
    addressing information.

SourceFiles
    ijkMesh.C
    ijkMeshI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ijkMesh_H
#define Foam_ijkMesh_H

#include "ijkAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class hexCell;

/*---------------------------------------------------------------------------*\
                           Class ijkMesh Declaration
\*---------------------------------------------------------------------------*/

class ijkMesh
:
    public ijkAddressing
{
public:

    // Constructors

        //- Construct zero-sized
        inline ijkMesh();

        //- Construct with addressing
        inline explicit ijkMesh(const labelVector& ijk);

        //- Construct with addressing
        inline ijkMesh(const label nx, const label ny, const label nz);


    // Member Functions

    // Access

        //- The number of mesh points (nx+1)*(ny+1)*(nz+1) in the i-j-k mesh.
        inline label nPoints() const;

        //- The number of mesh cells (nx*ny*nz) in the i-j-k mesh
        //  This is the same as the ijkAddressing::size()
        inline label nCells() const;

        //- The total number of mesh faces in the i-j-k mesh
        inline label nFaces() const;

        //- The number of internal faces in the i-j-k mesh
        inline label nInternalFaces() const;

        //- The number of boundary faces in the i-j-k mesh
        inline label nBoundaryFaces() const;

        //- The number of boundary faces on the box face
        //  Face: 0=x-min, 1=x-max, 2=y-min, 3=y-max, 4=z-min, 5=z-max
        inline label nBoundaryFaces(const direction shapeFacei) const;

        //- The linear cell index for an i-j-k position - same as index()
        inline label cellLabel
        (
            const label i,
            const label j,
            const label k
        ) const;

        //- The linear cell index for an i-j-k position - same as index()
        inline label cellLabel(const labelVector& ijk) const;

        //- The linear point index for an i-j-k position.
        //  Addressable in the range (ni+1, nj+1, nk+1).
        inline label pointLabel
        (
            const label i,
            const label j,
            const label k
        ) const;

        //- The linear point index for an i-j-k position.
        //  Addressable in the range (ni+1, nj+1, nk+1).
        inline label pointLabel(const labelVector& ijk) const;


    // Mesh Components

        //- The hex cell vertices for an i-j-k position
        hexCell vertLabels(const label i, const label j, const label k) const;

        //- The hex cell vertices for an i-j-k position
        hexCell vertLabels(const labelVector& ijk) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ijkMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
