/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::isoSurfaceCell

Description
    A surface formed by the iso value.
    After "Polygonising A Scalar Field Using Tetrahedrons", Paul Bourke
    (http://paulbourke.net/geometry/polygonise) and
    "Regularised Marching Tetrahedra: improved iso-surface extraction",
    G.M. Treece, R.W. Prager and A.H. Gee.

    See isoSurface. This is a variant which does tetrahedrisation from
    triangulation of face to cell centre instead of edge of face to two
    neighbouring cell centres. This gives much lower quality triangles
    but they are local to a cell.

SourceFiles
    isoSurfaceCell.C
    isoSurfaceCellTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceCell_H
#define isoSurfaceCell_H

#include "labelPair.H"
#include "pointIndexHit.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class triSurface;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceCell Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceCell
:
    public isoSurfaceBase
{
    // Private Data

        //- When to merge points
        const scalar mergeDistance_;

        //- The cell cut type
        List<cutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- For every unmerged triangle point the point in the triSurface
        labelList triPointMergeMap_;

        //- triSurface points that have weighted interpolation
        DynamicList<label> interpolatedPoints_;

        //- corresponding original, unmerged points
        DynamicList<FixedList<label, 3>> interpolatedOldPoints_;

        //- corresponding weights
        DynamicList<FixedList<scalar, 3>> interpolationWeights_;


    // Private Member Functions

        //- Get location of iso value as fraction inbetween s0,s1
        scalar isoFraction(const scalar s0, const scalar s1) const;

        //- Return the two common points between two triangles
        static labelPair findCommonPoints
        (
            const labelledTri&,
            const labelledTri&
        );

        //- Calculate centre of surface.
        static point calcCentre(const triSurface&);

        //- Replace surface (localPoints, localTris) with single point.
        //  Returns point.  Destroys arguments.
        pointIndexHit collapseSurface
        (
            const label celli,
            pointField& localPoints,
            DynamicList<labelledTri, 64>& localTris
        ) const;

        //- Determine per cc whether all near cuts can be snapped to single
        //- point.
        void calcSnappedCc
        (
            const bitSet& isTet,
            const scalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& triPoints,
            labelList& snappedCc
        ) const;

        //- Generate triangles for face connected to pointi
        void genPointTris
        (
            const scalarField& cellValues,
            const scalarField& pointValues,
            const label pointi,
            const label facei,
            const label celli,
            DynamicList<point, 64>& localTriPoints
        ) const;

        //- Generate triangles for tet connected to pointi
        void genPointTris
        (
            const scalarField& pointValues,
            const label pointi,
            const label facei,
            const label celli,
            DynamicList<point, 64>& localTriPoints
        ) const;

        //- Determine per point whether all near cuts can be snapped to single
        //  point.
        void calcSnappedPoint
        (
            const bitSet& isTet,
            const scalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& triPoints,
            labelList& snappedPoint
        ) const;

        //- Generate single point by interpolation or snapping
        template<class Type>
        Type generatePoint
        (
            const DynamicList<Type>& snappedPoints,
            const scalar s0,
            const Type& p0,
            const label p0Index,
            const scalar s1,
            const Type& p1,
            const label p1Index
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const DynamicList<Type>& snapped,
            const scalar s0,
            const Type& p0,
            const label p0Index,
            const scalar s1,
            const Type& p1,
            const label p1Index,
            const scalar s2,
            const Type& p2,
            const label p2Index,
            const scalar s3,
            const Type& p3,
            const label p3Index,
            DynamicList<Type>& points
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const scalarField& cVals,
            const scalarField& pVals,

            const Field<Type>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        triSurface stitchTriPoints
        (
            const bool checkDuplicates,
            const List<point>& triPoints,
            labelList& triPointReverseMap,  // unmerged to merged point
            labelList& triMap               // merged to unmerged triangle
        ) const;

        //- Determine edge-face addressing
        void calcAddressing
        (
            const triSurface& surf,
            List<FixedList<label, 3>>& faceEdges,
            labelList& edgeFace0,
            labelList& edgeFace1,
            Map<labelList>& edgeFacesRest
        ) const;

        //- Is triangle (given by 3 edges) not fully connected?
        static bool danglingTriangle
        (
            const FixedList<label, 3>& fEdges,
            const labelList& edgeFace1
        );

        //- Mark all non-fully connected triangles
        static label markDanglingTriangles
        (
            const List<FixedList<label, 3>>& faceEdges,
            const labelList& edgeFace0,
            const labelList& edgeFace1,
            const Map<labelList>& edgeFacesRest,
            boolList& keepTriangles
        );

        static triSurface subsetMesh
        (
            const triSurface&,
            const labelList& newToOldFaces,
            labelList& oldToNewPoints,
            labelList& newToOldPoints
        );


    // Sampling

        //- Interpolates cCoords, pCoords.
        template<class Type>
        tmp<Field<Type>> interpolateTemplate
        (
            const Field<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;

public:

    //- Runtime type information
    TypeName("isoSurfaceCell");


    // Constructors

        //- Construct from cell and point values
        //
        //  \param ignoreCells cells to ignore in the cellValues
        //
        //  Control parameters include
        //  - bounds optional bounding box for trimming
        //  - mergeTol fraction of mesh bounding box for merging points
        isoSurfaceCell
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bitSet& ignoreCells = bitSet()
        );


    //- Destructor
    virtual ~isoSurfaceCell() = default;


    // Sampling

        declareIsoSurfaceInterpolateMethod(scalar);
        declareIsoSurfaceInterpolateMethod(vector);
        declareIsoSurfaceInterpolateMethod(sphericalTensor);
        declareIsoSurfaceInterpolateMethod(symmTensor);
        declareIsoSurfaceInterpolateMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceCellTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
