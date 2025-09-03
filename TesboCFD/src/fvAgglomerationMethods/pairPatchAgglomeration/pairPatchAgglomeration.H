/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020,2023 OpenCFD Ltd.
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
    Foam::pairPatchAgglomeration

Description
    Primitive patch pair agglomerate method.

SourceFiles
    pairPatchAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pairPatchAgglomeration_H
#define Foam_pairPatchAgglomeration_H

#include "mathematicalConstants.H"
#include "indirectPrimitivePatch.H"
#include "List.H"
#include "edgeHashes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pairPatchAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class pairPatchAgglomeration
{
protected:

    typedef PrimitivePatch<List<face>, const pointField> bPatch;

    // Protected data

        //- Number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
        label mergeLevels_;

        //- Max number of levels
        label maxLevels_;

        //- Number of faces in coarsest level
        label nFacesInCoarsestLevel_;

        //- Global number of faces in coarsest level
        label nGlobalFacesInCoarsestLevel_;

        //- Feature angle
        scalar featureAngle_;

        //- The number of faces in each level
        labelList nFaces_;

        //- Cell restriction addressing array.
        //  Maps from the finer to coarse level
        PtrList<labelField> restrictAddressing_;

        //- Maps from finest to coarsest
        labelList restrictTopBottomAddressing_;

        //- Hierarchy of patch addressing
        PtrList<bPatch> patchLevels_;

        //- Edge weights
        EdgeMap<scalar> facePairWeight_;


private:

    // Private Member Functions

        //- Assemble coarse patch
        bool agglomeratePatch
        (
            const bPatch& patch,
            const labelList& fineToCoarse,
            const label fineLevelIndex,
            label& nMarkedEdges
        );

        //- Does combining facei with faceIDs produce a single face?
        bool isSingleEdgeLoop
        (
            const bPatch& patch,
            const labelList& faceIDs,
            const label facei
        ) const;

        //- Select neighbour with highest inbetween edge weight. Either looks
        //- at already clustered faces (addToCluster=true) or only
        //  unclustered
        label maxValidNeighbour
        (
            const bool addToCluster,
            const bPatch& patch,
            const label facei,
            const labelList& coarseCellMap
            //const labelListList& coarseToFine
        ) const;

        //- Agglomerate one level
        tmp<labelField> agglomerateOneLevel
        (
            label& nCoarseCells,
            const bPatch& patchLevel
        );

        //- Combine levels
        void combineLevels(const label curLevel);

        //- Shrink the number of levels to that specified
        void compactLevels(const label fineLevelIndex);

        //- Check the need for further agglomeration
        bool continueAgglomerating
        (
            const label nLocal,
            const label nLocalOld,
            const label nMarkedEdges
        );

        //- Set edge weights
        void setEdgeWeights(const label indexLevel);

        //- Set base patch edge weights
        void setLevel0EdgeWeights();

        //- Maps current level with base patch
        void mapBaseToTopAgglom(const label fineLevelIndex);

        //- No copy construct
        pairPatchAgglomeration(const pairPatchAgglomeration&) = delete;

        //- No copy assignment
        void operator=(const pairPatchAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("pairPatch");


    // Constructors

        //- Construct given faces, points and control dictionary
        pairPatchAgglomeration
        (
            const faceList& faces,
            const pointField& points,
            const dictionary& controlDict
        );

        //- Construct from components
        pairPatchAgglomeration
        (
            const faceList& faces,
            const pointField& points,
            const label mergeLevels,
            const label maxLevels,
            const label nFacesInCoarsestLevel,
            const label nGlobalFacesInCoarsestLevel,
            const scalar featureAngle = 0
        );


    // Destructor
    virtual ~pairPatchAgglomeration();


    // Member Functions

        //- Agglomerate patch
        void agglomerate();


        // Access

            //- Return size
            label size() const
            {
                return patchLevels_.size();
            }

            //- Return restriction from top level to bottom level
            const labelList& restrictTopBottomAddressing() const
            {
                return restrictTopBottomAddressing_;
            }

            //- Return primitivePatch of given level
            const bPatch& patchLevel(const label leveli) const;

            //- Return cell restrict addressing of given level
            const labelField& restrictAddressing(const label leveli) const
            {
                return restrictAddressing_[leveli];
            }


        // Restriction and prolongation

            //- Restrict (integrate by summation) cell field
            template<class Type>
            void restrictField
            (
                Field<Type>& cf,
                const Field<Type>& ff,
                const label fineLevelIndex
            ) const;

            //- Prolong (interpolate by injection) cell field
            template<class Type>
            void prolongField
            (
                Field<Type>& ff,
                const Field<Type>& cf,
                const label coarseLevelIndex
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pairPatchAgglomerationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
