/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fvMeshSubset

Description
    Holds a reference to the original mesh (the baseMesh)
    and optionally to a subset of that mesh (the subMesh)
    with mapping lists for points, faces, and cells.

    Can be constructed or reset to subset on the list of selected cells,
    which it creates as subMesh consisting only of the desired cells,
    with the mapping list for points, faces, and cells.

    Places all exposed internal faces into either
    - a user supplied patch
    - a newly created patch "oldInternalFaces"

    - reset() does coupled patch subsetting as well.
      If it detects a face on a coupled patch 'losing' its neighbour
      it will move the face into the oldInternalFaces patch.

    - if a user supplied patch is used, it is up to the destination
      patchField to handle exposed internal faces (mapping from face -1).
      If not provided the default is to assign the internalField. All the
      basic patch field types (e.g. fixedValue) will give a warning and
      preferably derived patch field types should be used that know how to
      handle exposed faces (e.g. use uniformFixedValue instead of fixedValue)

SourceFiles
    fvMeshSubset.C
    fvMeshSubsetI.H
    fvMeshSubsetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMeshSubset_H
#define Foam_fvMeshSubset_H

#include "fvMesh.H"
#include "pointMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "bitSet.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshSubset Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubset
{
    // Private Data

        //- The base mesh to subset from
        const fvMesh& baseMesh_;

        //- Demand-driven subset mesh (pointer)
        autoPtr<fvMesh> subMeshPtr_;

        //- Optional face mapping array with flip encoded (-1/+1)
        mutable autoPtr<labelList> faceFlipMapPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;

        //- PointPatch mapping array
        labelList pointPatchMap_;


    // Private Member Functions

        //- Modify nCellsUsingFace for coupled faces becoming 'uncoupled.
        void doCoupledPatches
        (
            const bool syncPar,
            labelUList& nCellsUsingFace
        ) const;

        //- Create zones for subMesh
        void subsetZones();

        //- Calculate face flip map
        void calcFaceFlipMap() const;


protected:

    // Protected Member Functions

        //- FatalError if subset has not been performed
        bool checkHasSubMesh() const;

        //- No copy construct
        fvMeshSubset(const fvMeshSubset&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubset&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("fvMeshSubset");


    // Static Data Members

        //- Name for exposed internal faces (default: oldInternalFaces)
        static word exposedPatchName;


    // Constructors

        //- Construct using the entire mesh (no subset)
        explicit fvMeshSubset(const fvMesh& baseMesh);

        //- Construct a zero-sized subset mesh, non-processor patches only
        fvMeshSubset(const fvMesh& baseMesh, const Foam::zero);

        //- Construct for a cell-subset of the given mesh
        //  See reset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const bitSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See reset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const labelUList& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See reset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const labelHashSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See reset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const label regioni,
            const labelUList& regions,
            const label patchID = -1,
            const bool syncPar = true
        );


    // Member Functions

    // Access

        //- Original mesh
        inline const fvMesh& baseMesh() const noexcept;

        //- Return baseMesh or subMesh, depending on the current state.
        inline const fvMesh& mesh() const noexcept;

        //- Have subMesh?
        inline bool hasSubMesh() const noexcept;

        //- Return reference to subset mesh
        inline const fvMesh& subMesh() const;

        //- Return reference to subset mesh
        inline fvMesh& subMesh();

        //- Return point map
        inline const labelList& pointMap() const;

        //- Return face map
        inline const labelList& faceMap() const;

        //- Return face map with sign to encode flipped faces
        inline const labelList& faceFlipMap() const;

        //- Return cell map
        inline const labelList& cellMap() const;

        //- Return patch map
        inline const labelList& patchMap() const;

        //- Return point-patch map. Usually identical to patchMap except if
        //- additional patches are added to the pointMesh.
        inline const labelList& pointPatchMap() const;


    // Edit

        //- Reset subMesh and all maps
        void clear();

        //- Reset subMesh and all maps. Same as clear()
        void reset();

        //- Reset to a zero-sized subset mesh, non-processor patches only
        void reset(const Foam::zero);

        //- Reset from components
        void reset
        (
            autoPtr<fvMesh>&& subMeshPtr,   //!< Mesh subset
            labelList&& pointMap,           //!< Point mapping
            labelList&& faceMap,            //!< Face mapping
            labelList&& cellMap,            //!< Cell mapping
            labelList&& patchMap            //!< Patch mapping
        );

        //- Use the specified subset of cells.
        //
        //  \par selectedCells The subset of cells to use
        //  \par patchID Patch id for exposed internal faces.
        //       Uses existing or creates "oldInternalFaces" for patchID == -1.
        //  \par syncPar
        //
        //  \note Handles coupled patches if necessary by making
        //      coupled patch faces part of patchID (ie, uncoupled)
        void reset
        (
            const bitSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Use the specified subset of cells.
        void reset
        (
            const labelUList& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Use the specified subset of cells.
        void reset
        (
            const labelHashSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Use the cells of cells corresponding to where region == regioni.
        void reset
        (
            const label regioni,
            const labelUList& regions,
            const label patchID = -1,
            const bool syncCouples = true
        );


    // Legacy method names (v2112 and earlier)

        //- Use the specified subset of cells. Same as reset()
        void setCellSubset
        (
            const bitSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        )
        {
            reset(selectedCells, patchID, syncPar);
        }

        //- Use the specified subset of cells. Same as reset()
        void setCellSubset
        (
            const labelUList& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        )
        {
            reset(selectedCells, patchID, syncPar);
        }

        //- Use the specified subset of cells. Same as reset()
        void setCellSubset
        (
            const labelHashSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        )
        {
            reset(selectedCells, patchID, syncPar);
        }

        //- Use the specified subset of cells. Same as reset()
        void setCellSubset
        (
            const label regioni,
            const labelUList& regions,
            const label patchID = -1,
            const bool syncPar = true
        )
        {
            reset(regioni, regions, patchID, syncPar);
        }


    // Field Mapping (static functions)

        //- Map volume internal (dimensioned) field
        template<class Type>
        static tmp<DimensionedField<Type, volMesh>>
        interpolate
        (
            const DimensionedField<Type, volMesh>&,
            const fvMesh& sMesh,
            const labelUList& cellMap
        );

        //- Map volume field.
        //  Optionally allow unmapped faces not to produce a warning
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const fvMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& cellMap,
            const labelUList& faceMap,
            const bool allowUnmapped = false
        );

        //- Map surface field.
        //  Optionally negates value if flipping a face
        //  (from exposing an internal face)
        template<class Type>
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const fvMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& cellMap,
            const labelUList& faceMap
        );

        //- Map point field
        template<class Type>
        static tmp<GeometricField<Type, pointPatchField, pointMesh>>
        interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            const pointMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& pointMap
        );


    // Field Mapping

        //- Map volume internal (dimensioned) field
        //- Optional unmapped argument (currently unused)
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        interpolate
        (
            const DimensionedField<Type, volMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map volume field.
        //  Optionally allow unmapped faces not to produce a warning
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map surface field.
        //  Optionally allow unmapped faces not to produce a warning
        //  (currently unused)
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map point field.
        //  Optionally allow unmapped points not to produce a warning
        //  (currently unused)
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            const bool allowUnmapped = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMeshSubsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
