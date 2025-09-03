/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::VF::voxel

Description
    Ray search engine based on discretising space into uniform voxels

    Voxels are refined using 2x2x2 refinement.

    The number of rays per face is supplied by the user, whereby rays are
    issued uniformly across a hemisphere.

Usage
    Minimal example by using \c <constant>/viewFactorsDict:
    \verbatim
    // Mandatory entries
    nRayPerFace        <label>;

    // Optional entries
    nTriPerVoxelMax     <label>;
    depthMax            <label>;

    // Inherited entries
    ...
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type | Reqd | Deflt
      nRayPerFace | Number of rays issued per face    | label | yes | -
      nRayPerFace | Target number of triangles per voxel | label | no | 50
      depthMax    | Maximum voxel refinement depth    | label | no  | 5
    \endtable

    The inherited entries are elaborated in:
    - \link raySearchEngine.H \endlink

SourceFiles
    voxelRaySearchEngine.C

SeeAlso
- Foam::VF::raySearchEngine

\*---------------------------------------------------------------------------*/

#ifndef Foam_vf_voxelRaySearchEngine_H
#define Foam_vf_voxelRaySearchEngine_H

#include "raySearchEngine.H"
#include "labelVector.H"
#include "OBJstream.H"
#include "pointIndexHit.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace VF
{

/*---------------------------------------------------------------------------*\
                            Class voxel Declaration
\*---------------------------------------------------------------------------*/

class voxel
:
    public raySearchEngine
{
    // Private Data

        //- Triangulation of view factor patches
        triSurface surface_;

        //- Triangle to mesh face addressing
        labelList triToGlobalFace_;

        //- Surface bounding box
        boundBox bb0_;

        //- Span of surface bounding box
        vector span0_;

        //- Voxel discretisation
        labelVector nijk_;

        //- Voxel dimensions
        vector dxyz_;

        //- Number of rays issued per face
        const label nRayPerFace_;

        //- Maximum number of triangles per voxel (trigger to refine voxels)
        const label nTriPerVoxelMax_;

        //- Maximum depth of voxel refinement. Note: voxels are refined 2x2x2
        const label depthMax_;

        //- List of triangle per voxel
        List<DynamicList<label>> objects_;

        //- List of triangle bounding boxes
        List<boundBox> objectBbs_;


    // Private Member Functions

        inline bool outOfBounds(const labelVector& ijk, const label dir) const;

        inline point localPosition(const vector& globalPosition) const;

        inline point globalPosition(const vector& localPosition) const;

        inline void setVoxelDims(const label i, const label j, const label k);

        inline void refineVoxelDims();

        inline point voxelMin
        (
            const label i,
            const label j,
            const label k
        ) const;

        inline point voxelMax
        (
            const label i,
            const label j,
            const label k
        ) const;

        inline constexpr label sign0(const scalar x) const;

        //- Set triangulation based on original mesh
        void setTriangulation(const fvMesh& mesh);

        //- Set the participating face vertices when simplifying edges
        static labelList setFaceVertexHits
        (
            const fvMesh& mesh,
            const labelList& patchIDs
        );

        //- Set triangulation based on agglomeration
        void setCoarseTriangulation(const fvMesh& mesh);

        //- Broadcast triangulation across all procs
        void broadcast();

        void refineObjects
        (
            List<DynamicList<label>>& objects,
            const label triMax
        );

        label addBbToVoxels
        (
            const boundBox& bb,
            const label trii,
            List<DynamicList<label>>& objects
        ) const;

        void voxelise
        (
            List<DynamicList<label>>& objects,
            const label trii0,
            const label depth
        );

        pointHit rayTriIntersect
        (
            const label trii,
            const point& origin,
            const vector& direction
        ) const;

        pointIndexHit hitObject
        (
            const label voxeli,
            const point& origin,
            const vector& direction,
            const vector& t,
            const scalar minDistance = 1e-6
        ) const;

        void writeBox
        (
            OBJstream& os,
            bool lines,
            const boundBox& bb
        ) const;

        void writeVoxels(const word& fName) const;

        void writeTriBoundBoxes(const word& fName) const;

        void writeHitObjects
        (
            const label voxeli,
            const point& origin,
            const vector& dir
        ) const;


public:

    //- Runtime type information
    TypeName("voxel");

    //- Constructor
    voxel(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~voxel() = default;


    // Public Member Functions

        inline labelVector nijk() const noexcept;

        inline label nVoxel() const noexcept;

        inline label voxeli(const labelVector ijk) const noexcept;

        inline label voxeli
        (
            const label i,
            const label j,
            const label k
        ) const noexcept;

        inline labelVector ijk(const label voxeli) const noexcept;

        pointIndexHit hit(const point& origin, const vector& dir) const;

        pointIndexHit hit(const label tri0, const vector& dir) const;

        virtual void shootRays
        (
            labelList& rayStartFaceOut,
            labelList& rayEndFaceOut
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "voxelRaySearchEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
