/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::snappyVoxelMeshDriver

Description
    Equivalent of snappyRefineDriver but operating on a voxel mesh.

    Used to estimate cell size count from refinement (currently).

SourceFiles
    snappyVoxelMeshDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyVoxelMeshDriver_H
#define snappyVoxelMeshDriver_H

#include "labelList.H"
#include "DynamicList.H"
#include "labelVector.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class refinementParameters;
class layerParameters;
class meshRefinement;

/*---------------------------------------------------------------------------*\
                           Class snappyVoxelMeshDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyVoxelMeshDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;

        //- Bounding box
        boundBox bb_;

        //- Number of voxels
        labelVector n_;

        //- Per voxel the refinement level
        labelList voxelLevel_;

        //- Per insersected voxel the originating surface
        labelList globalRegion_;


    // Private Member Functions

        void addNeighbours
        (
            const labelList& cellLevel,
            const labelVector& voxel,
            const label voxeli,
            DynamicList<labelVector>& front
        ) const;

        tmp<pointField> voxelCentres() const;

        void isInside
        (
            const pointField& voxelCentres,
            boolList& isVoxelInMesh
        ) const;

        void markSurfaceRefinement
        (
            labelList& voxelLevel,
            labelList& globalRegion
        ) const;

        void findVoxels
        (
            const labelList& voxelLevel,
            const pointField& locationsOutsideMesh,
            labelList& voxels
        ) const;

        void floodFill
        (
            const label voxeli,
            const label newLevel,
            labelList& voxelLevel
        ) const;

        void max
        (
            const labelList& maxLevel,
            labelList& voxelLevel
        ) const;

        labelList count(const labelList& voxelLevel) const;


       //- No copy construct
        snappyVoxelMeshDriver(const snappyVoxelMeshDriver&) = delete;

        //- No copy assignment
        void operator=(const snappyVoxelMeshDriver&) = delete;


public:

    //- Runtime type information
    ClassName("snappyVoxelMeshDriver");


    // Constructors

        //- Construct from components
        snappyVoxelMeshDriver
        (
            meshRefinement& meshRefiner,
            const labelUList& globalToMasterPatch,
            const labelUList& globalToSlavePatch
        );


    // Member Functions

        void doRefine(const refinementParameters& refineParams);

        //void doLayers(const layerParameters& layerParams);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
