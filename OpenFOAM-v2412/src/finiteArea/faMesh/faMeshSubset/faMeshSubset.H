/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::faMeshSubset

Description
    Holds a reference to the original mesh (the baseMesh)
    and optionally to a subset of that mesh (the subMesh)
    with mapping lists for points, faces, and cells.

Caution
    Currently not really functional for subsetting beyond handling
    a simple zero-sized mesh.

SourceFiles
    faMeshSubset.C
    faMeshSubsetI.H
    faMeshSubsetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshSubset_H
#define Foam_faMeshSubset_H

#include "areaFaMesh.H"
#include "edgeFaMesh.H"
#include "GeometricField.H"
#include "bitSet.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                        Class faMeshSubset Declaration
\*---------------------------------------------------------------------------*/

class faMeshSubset
{
    // Private Data

        //- The base mesh to subset from
        const faMesh& baseMesh_;

        //- Demand-driven subset mesh (pointer)
        autoPtr<faMesh> subMeshPtr_;

        //- Optional edge mapping array with flip encoded (-1/+1)
        mutable autoPtr<labelList> edgeFlipMapPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;


    // Private Member Functions

        //- Calculate face flip map
        void calcEdgeFlipMap() const;


protected:

    // Protected Member Functions

        //- FatalError if subset has not been performed
        bool checkHasSubMesh() const;


public:

    // Static Data Members

        //- Name for exposed internal edges (default: oldInternalEdges)
        static word exposedPatchName;

    // Generated Methods

        //- No copy construct
        faMeshSubset(const faMeshSubset&) = delete;

        //- No copy assignment
        void operator=(const faMeshSubset&) = delete;


    // Constructors

        //- Construct using the entire mesh (no subset)
        explicit faMeshSubset(const faMesh& baseMesh);

        //- Construct a zero-sized subset mesh, non-processor patches only
        faMeshSubset(const faMesh& baseMesh, const Foam::zero);


    // Member Functions

    // Access

        //- Original mesh
        inline const faMesh& baseMesh() const noexcept;

        //- Return baseMesh or subMesh, depending on the current state.
        inline const faMesh& mesh() const noexcept;

        //- Have subMesh?
        inline bool hasSubMesh() const noexcept;

        //- Return reference to subset mesh
        inline const faMesh& subMesh() const;

        //- Return reference to subset mesh
        inline faMesh& subMesh();

        //- Return point map
        inline const labelList& pointMap() const;

        //- Return face map
        inline const labelList& faceMap() const;

        //- Return edge map with sign to encode flipped edges
        inline const labelList& edgeFlipMap() const;

        //- Return cell map
        inline const labelList& cellMap() const;

        //- Return patch map
        inline const labelList& patchMap() const;


    // Edit

        //- Reset subMesh and all maps
        void clear();

        //- Reset subMesh and all maps. Same as clear()
        void reset();

        //- Reset to a zero-sized subset mesh, non-processor patches only
        void reset(const Foam::zero);


    // Field Mapping (static functions)

        //- Map area field.
        //  Optionally allow unmapped faces not to produce a warning
        template<class Type>
        static tmp<GeometricField<Type, faPatchField, areaMesh>>
        interpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const faMesh& sMesh,
            const bool allowUnmapped = false
        );

        //- Map edge field.
        //  Optionally allow unmapped faces not to produce a warning
        template<class Type>
        static tmp<GeometricField<Type, faePatchField, edgeMesh>>
        interpolate
        (
            const GeometricField<Type, faePatchField, edgeMesh>&,
            const faMesh& sMesh
        );


    // Field Mapping

        //- Map area field.
        //  Optionally allow unmapped faces not to produce a warning
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        interpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map edge field.
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        interpolate
        (
            const GeometricField<Type, faePatchField, edgeMesh>&,
            const bool allowUnmapped = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faMeshSubsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faMeshSubsetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
