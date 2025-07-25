/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::cuttingPlane

Description
    Constructs cutting plane through a mesh.

    No attempt at resolving degenerate cases.
    Since the cut faces can be quite ugly, they will often be triangulated.

SourceFiles
    cuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPlane_H
#define cuttingPlane_H

#include "plane.H"
#include "cuttingSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class cuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlane
:
    public plane,
    public cuttingSurfaceBase
{
    // Private Member Functions

        //- Determine cut cells, possibly restricted to a list of cells
        //
        //  \param sides [out] For each mesh point, the encoded side of the
        //      plane (0=BACK, 1=ONPLANE, 2=FRONT).
        //  \param cellCuts [in,out] On input an empty set (ie, no restriction)
        //      or subsetted cells. On output, the cells cut according to the
        //      plane-sides detection.
        //
        //  \return number of faces cut
        label calcCellCuts
        (
            const primitiveMesh& mesh,
            PackedList<2>& sides,
            bitSet& cellCuts
        );


protected:

    // Constructors

        //- Construct from a plane description without any cutting
        cuttingPlane(const plane& pln);


    // Protected Member Functions

        //- Cut mesh, restricted to a list of cells
        using cuttingSurfaceBase::performCut;

        //- Cut mesh, restricted to a list of cells
        //  Reclaim memory for cellIdLabels
        virtual void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellIdLabels
        );


        //- Check and warn if bounding boxes do not intersect,
        //- and if the plane does not intersect the bounding boxes
        void checkOverlap
        (
            const word callerName,
            const boundBox& meshBounds,
            const boundBox& userBounds
        ) const;

        //- Define cell selection from bounding-box and zones.
        //
        //  \param userBounds Optionally user-specified bounding box
        //  \param zoneNames  Optionally user-specified zone names
        //  \param warn Check and warn if the plane does not
        //      intersect with the bounds of the mesh (or submesh) or
        //      if the bounding box does not overlap with the mesh (or submesh)
        //
        //  \return A set of nCells size with the selected cells or an empty
        //      set if no bounding-box or zones were specified.
        bitSet cellSelection
        (
            const polyMesh& mesh,
            const boundBox& userBounds,
            const wordRes& zoneNames,
            const word callerName,
            const bool warn
        ) const;


public:

    // Constructors

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            const bitSet& cellIdLabels = bitSet()
        );

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellIdLabels
        );

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& cellIdLabels
        );


    // Member Functions

        //- The plane used.
        const plane& planeDesc() const
        {
            return *this;
        }

        //- The plane used. Non-const access.
        plane& planeDesc()
        {
            return *this;
        }


    // Member Operators

        //- Copy assignment
        void operator=(const cuttingPlane& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
