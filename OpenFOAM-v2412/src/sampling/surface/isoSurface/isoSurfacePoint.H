/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::isoSurfacePoint

Description
    A surface formed by the iso value.
    After "Regularised Marching Tetrahedra: improved iso-surface extraction",
    G.M. Treece, R.W. Prager and A.H. Gee.

    Note:
    - does tets without using cell centres/cell values. Not tested.
    - regularisation can create duplicate triangles/non-manifold surfaces.
    Current handling of those is bit ad-hoc for now and not perfect.
    - regularisation does not do boundary points so as to preserve the
      boundary perfectly.
    - uses geometric merge with fraction of bounding box as distance.
    - triangles can be between two cell centres so constant sampling
      does not make sense.
    - on empty patches behaves like zero gradient.
    - does not do 2D correctly, creates non-flat iso surface.
    - on processor boundaries might two overlapping (identical) triangles
      (one from either side)

    The handling on coupled patches is a bit complex. All fields
    (values and coordinates) get rewritten so
    - empty patches get zerogradient (value) and facecentre (coordinates)
    - separated processor patch faces get interpolate (value) and facecentre
      (coordinates). (this is already the default for cyclics)
    - non-separated processor patch faces keep opposite value and cell centre

    Now the triangle generation on non-separated processor patch faces
    can use the neighbouring value. Any separated processor face or cyclic
    face gets handled just like any boundary face.

SourceFiles
    isoSurfacePoint.C
    isoSurfacePointTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_isoSurfacePoint_H
#define Foam_isoSurfacePoint_H

#include "bitSet.H"
#include "volFields.H"
#include "slicedVolFields.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class plane;
class treeBoundBox;
class triSurface;

/*---------------------------------------------------------------------------*\
                       Class isoSurfacePoint Declaration
\*---------------------------------------------------------------------------*/

class isoSurfacePoint
:
    public isoSurfaceBase
{
    // Private Data

        //- Cell values.
        //- Input volScalarField with separated coupled patches rewritten
        autoPtr<slicedVolScalarField> cValsPtr_;

        //- When to merge points
        const scalar mergeDistance_;

        //- Whether face might be cut
        List<cutType> faceCutType_;

        //- Whether cell might be cut
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

        // Point synchronisation

            //- Per face whether is collocated
            static bitSet collocatedFaces(const coupledPolyPatch& cpp);

            //- Synchronise points on all non-separated coupled patches
            void syncUnseparatedPoints
            (
                pointField& collapsedPoint,
                const point& nullValue
            ) const;


        //- Get location of iso value as fraction inbetween s0,s1
        scalar isoFraction
        (
            const scalar s0,
            const scalar s1
        ) const;

        //- Check if any edge of a face is cut
        bool isEdgeOfFaceCut
        (
            const scalarField& pVals,
            const face& f,
            const bool ownLower,
            const bool neiLower
        ) const;

        //- Get neighbour value and position.
        void getNeighbour
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const label celli,
            const label facei,
            scalar& nbrValue,
            point& nbrPoint
        ) const;

        //- Determine for every face/cell whether it (possibly) generates
        //  triangles.
        void calcCutTypes
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals
        );

        static point calcCentre(const triSurface&);

        //- Determine per cc whether all near cuts can be snapped to single
        //  point.
        void calcSnappedCc
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedCc
        ) const;

        //- Determine per point whether all near cuts can be snapped to single
        //  point.
        void calcSnappedPoint
        (
            const bitSet& isBoundaryPoint,
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedPoint
        ) const;


        //- Return input field with coupled (and empty) patch values rewritten
        template<class Type>
        tmp<VolumeSliceField<Type>>
        adaptPatchFields(const VolumeField<Type>& fld) const;

        //- Generate single point by interpolation or snapping
        template<class Type>
        Type generatePoint
        (
            const scalar s0,
            const Type& p0,
            const bool hasSnap0,
            const Type& snapP0,

            const scalar s1,
            const Type& p1,
            const bool hasSnap1,
            const Type& snapP1
        ) const;


        //- Note: cannot use simpler isoSurfaceCell::generateTriPoints since
        //  the need here to sometimes pass in remote 'snappedPoints'
        template<class Type>
        void generateTriPoints
        (
            const scalar s0,
            const Type& p0,
            const bool hasSnap0,
            const Type& snapP0,

            const scalar s1,
            const Type& p1,
            const bool hasSnap1,
            const Type& snapP1,

            const scalar s2,
            const Type& p2,
            const bool hasSnap2,
            const Type& snapP2,

            const scalar s3,
            const Type& p3,
            const bool hasSnap3,
            const Type& snapP3,

            DynamicList<Type>& pts
        ) const;

        template<class Type>
        label generateFaceTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const VolumeField<Type>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,
            const label facei,

            const scalar neiVal,
            const Type& neiPt,
            const bool hasNeiSnap,
            const Type& neiSnapPt,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const VolumeField<Type>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        template<class Type>
        static tmp<Field<Type>> interpolate
        (
            const label nPoints,
            const labelList& triPointMergeMap,
            const labelList& interpolatedPoints,
            const List<FixedList<label, 3>>& interpolatedOldPoints,
            const List<FixedList<scalar, 3>>& interpolationWeights,
            const DynamicList<Type>& unmergedValues
        );

        triSurface stitchTriPoints
        (
            const bool checkDuplicates,
            const List<point>& triPoints,
            labelList& triPointReverseMap,  // unmerged to merged point
            labelList& triMap               // merged to unmerged triangle
        ) const;

        //- Trim triangle to planes
        static void trimToPlanes
        (
            const PtrList<plane>& planes,
            const triPointRef& tri,
            DynamicList<point>& newTriPoints
        );

        //- Trim all triangles to box
        static void trimToBox
        (
            const treeBoundBox& bb,
            DynamicList<point>& triPoints,
            DynamicList<label>& triMeshCells
        );

        //- Trim all triangles to box. Determine interpolation
        //  for existing and new points
        static void trimToBox
        (
            const treeBoundBox& bb,
            DynamicList<point>& triPoints,
            DynamicList<label>& triMap,
            labelList& triPointMap,
            labelList& interpolatedPoints,
            List<FixedList<label, 3>>& interpolatedOldPoints,
            List<FixedList<scalar, 3>>& interpolationWeights
        );

        static triSurface subsetMesh
        (
            const triSurface&,
            const labelList& newToOldFaces,
            labelList& oldToNewPoints,
            labelList& newToOldPoints
        );


        //- Interpolates cCoords, pCoords.
        //  Uses the references to the original fields used to create the
        //  iso surface.
        template<class Type>
        tmp<Field<Type>> interpolateTemplate
        (
            const VolumeField<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;


public:

    //- Declare friendship to share some functionality
    friend class isoSurfaceCell;
    friend class isoSurfaceTopo;


    //- Runtime type information
    TypeName("isoSurfacePoint");


    // Constructors

        //- Construct from cell values and point values.
        //  Uses boundaryField for boundary values.
        //  Holds reference to cellIsoVals and pointIsoVals.
        //
        //  Control parameters include
        //  - bounds optional bounding box for trimming
        //  - mergeTol fraction of mesh bounding box for merging points
        isoSurfacePoint
        (
            const volScalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bitSet& ignoreCells = bitSet()  //!< unused
        );


    //- Destructor
    virtual ~isoSurfacePoint() = default;


    // Member Functions

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
    #include "isoSurfacePointTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
