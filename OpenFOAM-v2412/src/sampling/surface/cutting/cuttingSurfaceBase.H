/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::cuttingSurfaceBase

Description
    Base for creating a MeshedSurface by performing some type of cell
    cutting/intersection.

    No attempt at resolving degenerate cases.
    Since the cut faces can be quite ugly, they will often be triangulated.

SourceFiles
    cuttingSurfaceBase.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingSurfaceBase_H
#define cuttingSurfaceBase_H

#include "bitSet.H"
#include "faceList.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                        Class cuttingSurfaceBase Declaration
\*---------------------------------------------------------------------------*/

class cuttingSurfaceBase
:
    public meshedSurface
{
protected:

    //- Typedef for convenience
    typedef meshedSurface Mesh;


    // Protected Data

        //- List of the cells cut
        labelList meshCells_;


    // Protected Member Functions

        //- Walk cell cuts to create faces
        //
        //  \tparam EdgeOrientIntersect
        //      Parameter (edge&), returns bool.
        //      Orient edge for a consistent positive gradient.
        //      Checks for edge intersection (true|false).
        //
        //  \tparam EdgeAlphaIntersect
        //      Parameter (const edge&), returns scalar.
        //      Determine alpha [0-1] for an intersecting edge.
        //      No guarantees when used with non-intersecting edges.
        //
        //  \param cellCuts [in] The cells to walk.
        template<class EdgeOrientIntersect, class EdgeAlphaIntersect>
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const bitSet& cellCuts,
            const EdgeOrientIntersect& edgeOrientIntersect,
            const EdgeAlphaIntersect&  edgeAlphaIntersect,
            const bool triangulate,
            label nFaceCuts = 0
        );


        //- Cut mesh, restricted to a list of cells
        virtual void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& cellIdLabels
        );

        //- Cut mesh, restricted to a list of cells
        virtual void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const bitSet& cellSelectionMask = bitSet()
        );

        //- Cut mesh, restricted to a list of cells
        //  Reclaim memory for cellSelectionMask
        virtual void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellSelectionMask
        ) = 0;

        //- Remap action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


        //- Check and warn if bounding boxes do not intersect
        static void checkOverlap
        (
            const word callerName,
            const boundBox& meshBounds,
            const boundBox& userBounds
        );


        //- Define cell selection from bounding-box and zones.
        //
        //  \param userBounds Optionally user-specified bounding box
        //  \param zoneNames  Optionally user-specified zone names
        //  \param meshBounds [out] The effective mesh bounds after applying
        //      the user-specified zone names
        //
        //  \return A set of nCells size with the selected cells or an empty
        //      set if no bounding-box or zones were specified.
        static bitSet cellSelection
        (
            const polyMesh& mesh,
            const boundBox& userBounds,
            const wordRes& zoneNames,
            boundBox& meshBounds
        );

        //- Define cell selection from bounding-box and zones.
        //
        //  \param userBounds Optionally user-specified bounding box
        //  \param zoneNames  Optionally user-specified zone names
        //  \param callerName The caller name for warnings
        //  \param warn Check and warn if the bounding box does not
        //      overlap with the mesh (or submesh)
        //
        //  \return A set of nCells size with the selected cells or an empty
        //      set if no bounding-box or zones were specified.
        static bitSet cellSelection
        (
            const polyMesh& mesh,
            const boundBox& userBounds,
            const wordRes& zoneNames,
            const word callerName,
            const bool warn
        );


public:

    //- Debug information
    static int debug;


    // Constructors

        //- Default construct
        cuttingSurfaceBase() = default;


    //- Destructors
    virtual ~cuttingSurfaceBase() = default;


    // Member Functions

        //- The mesh cells cut
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- The mesh cells cut
        labelList& meshCells()
        {
            return meshCells_;
        }

        //- Have any cells been cut?
        bool cut() const
        {
            return meshCells_.size();
        }


    // Member Operators

        //- Copy assignment
        void operator=(const cuttingSurfaceBase& rhs);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cuttingSurfaceBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
