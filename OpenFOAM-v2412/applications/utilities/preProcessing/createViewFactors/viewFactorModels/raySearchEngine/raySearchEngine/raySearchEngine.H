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
    Foam::VF::raySearchEngine

Description
    Base class for ray search engines

    Participating patches must be in the \c viewFactorWall group, i.e. using the
    \c inGroups entry of the "\<case\>/polyMesh/boundary" file.

    \verbatim
    myPatch
    {
        type            wall;
        inGroups        2(wall viewFactorWall);
        ...
    }
    \endverbatim

    Face agglomeration can be employed, created using the \c faceAgglomerate
    utility. The file name to be read can be user-defined:

    \verbatim
    // Name of agglomeration file; default = finalAgglom
    agglom      finalAgglom;
    \endverbatim

SourceFiles
    raySearchEngine.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vf_raySearchEngine_H
#define Foam_vf_raySearchEngine_H

#include "cartesianCS.H"
#include "mapDistribute.H"
#include "singleCellFvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace VF
{

/*---------------------------------------------------------------------------*\
                       Class raySearchEngine Declaration
\*---------------------------------------------------------------------------*/

class raySearchEngine
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Parallel map
        mutable autoPtr<mapDistribute> mapPtr_;

        //- Compact to global addressing
        mutable labelList compactToGlobal_;

        //- Global numbering
        globalIndex globalNumbering_;

        //- Name of patch group to identify participating patches
        const word patchGroup_;

        //- List of participating patch IDs
        labelList patchIDs_;

        //- Patch areas
        scalarList patchAreas_;

        //- Agglomeration flag
        bool agglomerate_;

        //- Agglomerated mesh representation
        autoPtr<singleCellFvMesh> agglomMeshPtr_;

        //- Number of original faces
        label nFace_;

        //- Number of coarse faces
        label nCoarseFace_;

        //- List of all face centres per processor
        List<pointField> allCf_;

        //- List of all face areas per processor
        List<vectorField> allSf_;

        //- List of all face agglomeration index per processor
        List<labelField> allAgg_;


    // Protected Member Functions

        static void check(const labelList& nVisibleFaceFaces);

        static label closestPointIndex
        (
            const point& p0,
            const List<point>& pts
        );

        //- Create patch geometry based on the original mesh
        void createGeometry();

        //- Create parallel addressing - map, compact-to-global
        void createParallelAddressing(labelList& rayEndFace) const;

        //- Create Cartesian co-ordinate system
        coordSystem::cartesian createCoordSystem
        (
            const point& origin,
            const vector& dir
        ) const;

        //- Create patch geometry based on the agglomerated mesh
        void createAgglomeration(const IOobject& io);

        //- Create a set of points describing a hemisphere
        //  Note: origin is (0 0 0)
        tmp<pointField> createHemiPoints(const label nRayPerFace) const;


public:

    static const label maxDynListLength;

    //- Run-time type information
    TypeName("raySearchEngine");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        raySearchEngine,
        mesh,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    //- Selector
    static autoPtr<raySearchEngine> New
    (
        const fvMesh& mesh,
        const dictionary& dict
    );


    // Generated Methods

        //- No copy construct
        raySearchEngine(const raySearchEngine&) = delete;

        //- No copy assignment
        void operator=(const raySearchEngine&) = delete;


    //- Constructor
    raySearchEngine(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~raySearchEngine() = default;


    // Public Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const noexcept;

            //- Parallel map
            inline const mapDistribute& map() const;

            //- Compact to global addressing
            inline const labelList& compactToGlobal() const noexcept;

            //- Global numbering
            inline const globalIndex& globalNumbering() const noexcept;

            //- List of participating patch IDs
            inline const labelList& patchIDs() const noexcept;

            //- Patch areas
            inline const scalarList& patchAreas() const noexcept;

            //- Number of participating faces
            inline label nParticipatingFaces() const;

            //- List of all face centres per processor
            inline const List<pointField>& allCf() const noexcept;

            //- List of all face areas per processor
            inline const List<vectorField>& allSf() const noexcept;

            //- List of all face agglomeration index per processor
            inline const List<labelField>& allAgg() const noexcept;


    // Main calculation functions

        //- Shoot rays; returns lists of ray start and end faces
        virtual void shootRays
        (
            labelList& rayStartFaceOut,
            labelList& rayEndFaceOut
        ) const = 0;

        //- Correct
        virtual void correct(labelListList& visibleFaceFaces) const;

        //- Create compact addressing
        void compactAddressing
        (
            const mapDistribute& map,
            pointField& compactCf,
            vectorField& compactSf,
            List<List<vector>>& compactFineSf,
            List<List<point>>& compactFineCf,
            DynamicList<List<point>>& compactPoints,
            DynamicList<label>& compactPatchId
        ) const;

        //- Interpolate field
        template<class Type>
        void interpolate
        (
            GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<Type>>& values
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "raySearchEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "raySearchEngineTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
