/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::AMIInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitive patches with an arbitrary mesh interface (AMI).

    Based on the algorithm given in:

        Conservative interpolation between volume meshes by local Galerkin
        projection, Farrell PE and Maddison JR, 2011, Comput. Methods Appl.
        Mech Engrg, Volume 200, Issues 1-4, pp 89-100

    Interpolation requires that the two patches should have opposite
    orientations (opposite normals).  The 'reverseTarget' flag can be used to
    reverse the orientation of the target patch.

SourceFiles
    AMIInterpolationI.H
    AMIInterpolation.C
    AMIInterpolationNew.C
    AMIInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_AMIInterpolation_H
#define Foam_AMIInterpolation_H

#include "searchableSurface.H"
#include "boolList.H"
#include "primitivePatch.H"
#include "ops.H"
#include "refPtr.H"
#include "pointList.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class AMIInterpolation Declaration
\*---------------------------------------------------------------------------*/

class AMIInterpolation
{
public:

    // Public Data Types

        static bool cacheIntersections_;


protected:

    //- Local typedef to octree tree-type
    typedef treeDataPrimitivePatch<primitivePatch> treeType;

    // Protected data

        //- Flag to indicate that the two patches must be matched/an overlap
        //- exists between them
        bool requireMatch_;

        //- Flag to indicate that the two patches are co-directional and
        //- that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Threshold weight below which interpolation is deactivated
        const scalar lowWeightCorrection_;

        //- Index of processor that holds all of both sides.
        //- The value is -1 for distributed cases
        label singlePatchProc_;

        //- Communicator to use for parallel operations.
        label comm_;


        // Source patch

            //- Source face areas
            scalarList srcMagSf_;

            //- Addresses of target faces per source face
            labelListList srcAddress_;

            //- Weights of target faces per source face
            scalarListList srcWeights_;

            //- Sum of weights of target faces per source face
            scalarField srcWeightsSum_;

            //- Centroid of target faces per source face
            pointListList srcCentroids_;

            //- Source patch points if input points are manipulated, e.g.
            //- projected
            pointField srcPatchPts_;

            //- Source patch using manipulated input points
            refPtr<primitivePatch> tsrcPatch0_;

            //- Source map pointer - parallel running only
            autoPtr<mapDistribute> srcMapPtr_;



        // Target patch

            //- Target face areas
            scalarList tgtMagSf_;

            //- Addresses of source faces per target face
            labelListList tgtAddress_;

            //- Weights of source faces per target face
            scalarListList tgtWeights_;

            //- Sum of weights of source faces per target face
            scalarField tgtWeightsSum_;

            //- Centroid of source faces per target face
            pointListList tgtCentroids_;

            //- Target patch points if input points are manipulated, e.g.
            //- projected
            pointField tgtPatchPts_;

            //- Target patch using manipulated input points
            refPtr<primitivePatch> ttgtPatch0_;

            //- Target map pointer - parallel running only
            autoPtr<mapDistribute> tgtMapPtr_;

        //- Up-to-date flag
        bool upToDate_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const AMIInterpolation&) = delete;


        // Initialisation

            //- Reset the octree for the patch face search
            autoPtr<indexedOctree<treeType>> createTree
            (
                const primitivePatch& patch
            ) const;

            //- Calculate if patches are on multiple processors
            label calcDistribution
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            ) const;

            //- Project points to surface
            void projectPointsToSurface
            (
                const searchableSurface& surf,
                pointField& pts
            ) const;


        // Access

            //- Return the orginal src patch with optionally updated points
            inline const primitivePatch& srcPatch0() const;

            //- Return the orginal tgt patch with optionally updated points
            inline const primitivePatch& tgtPatch0() const;


        // Evaluation

            //- Normalise the (area) weights - suppresses numerical error in
            //- weights calculation
            //  NOTE: if area weights are incorrect by 'a significant amount'
            //     normalisation may stabilise the solution, but will introduce
            //     numerical error!
            static void normaliseWeights
            (
                const scalarList& patchAreas,
                const word& patchName,
                const labelListList& addr,
                scalarListList& wght,
                scalarField& wghtSum,
                const bool conformal,
                const bool output,
                const scalar lowWeightTol,
                const label comm
            );


        // Constructor helpers

            static void agglomerate
            (
                const autoPtr<mapDistribute>& targetMap,
                const scalarList& fineSrcMagSf,
                const labelListList& fineSrcAddress,
                const scalarListList& fineSrcWeights,

                const labelList& sourceRestrictAddressing,
                const labelList& targetRestrictAddressing,

                scalarList& srcMagSf,
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarField& srcWeightsSum,
                autoPtr<mapDistribute>& tgtMap,
                const label comm
            );


public:

    //- Runtime type information
    TypeName("AMIInterpolation");

    // Selection tables

        //- Selection table for dictionary construction
        declareRunTimeSelectionTable
        (
            autoPtr,
            AMIInterpolation,
            dict,
            (
                const dictionary& dict,
                const bool reverseTarget
            ),
            (
                dict,
                reverseTarget
            )
        );

        //- Selection table for component-wise construction
        declareRunTimeSelectionTable
        (
            autoPtr,
            AMIInterpolation,
            component,
            (
                const bool requireMatch,
                const bool reverseTarget,
                const scalar lowWeightCorrection
            ),
            (
                requireMatch,
                reverseTarget,
                lowWeightCorrection
            )
        );

        //- Selector for dictionary
        static autoPtr<AMIInterpolation> New
        (
            const word& modelName,
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Selector for components
        static autoPtr<AMIInterpolation> New
        (
            const word& modelName,
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1
        );


    // Constructors

        //- Construct from dictionary
        explicit AMIInterpolation
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        explicit AMIInterpolation
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1
        );

        //- Construct from agglomeration of AMIInterpolation. Agglomeration
        //- passed in as new coarse size and addressing from fine from coarse
        AMIInterpolation
        (
            const AMIInterpolation& fineAMI,
            const labelList& sourceRestrictAddressing,
            const labelList& neighbourRestrictAddressing
        );

        //- Construct as copy
        AMIInterpolation(const AMIInterpolation& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>::New(*this);
        }

        //- Construct from Istream
        //  Note: does not transfer src(tgt)PatchPts, tsrc(tgt)Patch0_
        AMIInterpolation(Istream& is);


    //- Destructor
    virtual ~AMIInterpolation() = default;


    // Member Functions

        // Access

            //- Is up-to-date?
            bool upToDate() const noexcept { return upToDate_; }

            //- Set as up-to-date, return the previous value
            bool upToDate(bool flag) noexcept
            {
                bool old(upToDate_);
                upToDate_ = flag;
                return old;
            }

            //- Distributed across processors (singlePatchProc == -1)
            inline bool distributed() const noexcept;

            //- Return the requireMatch flag
            inline bool requireMatch() const noexcept;

            //- Set the require match flag, return the \em new value
            inline bool setRequireMatch(bool flag) noexcept;

            //- Return true if requireMatch and but not lowWeightCorrection
            inline bool mustMatchFaces() const;

            //- Access to the reverseTarget flag
            inline bool reverseTarget() const noexcept;

            //- Threshold weight below which interpolation is deactivated
            inline scalar lowWeightCorrection() const;

            //- Return true if employing a 'lowWeightCorrection'
            inline bool applyLowWeightCorrection() const;

            //- The processor holding all faces (both sides),
            //- or -1 if distributed
            inline label singlePatchProc() const noexcept;

            //- Communicator to use for parallel operations
            inline label comm() const;

            //- Set communicator to use for parallel operations. Return
            //  old value
            inline label comm(const label newComm);


            // Source patch

                //- Return const access to source patch face areas
                inline const List<scalar>& srcMagSf() const;

                //- Return access to source patch face areas
                inline List<scalar>& srcMagSf();

                //- Return const access to source patch addressing
                inline const labelListList& srcAddress() const;

                //- Return access to source patch addressing
                inline labelListList& srcAddress();

                //- Return const access to source patch weights
                inline const scalarListList& srcWeights() const;

                //- Return access to source patch weights
                inline scalarListList& srcWeights();

                //- Return const access to normalisation factor of source
                //- patch weights (i.e. the sum before normalisation)
                inline const scalarField& srcWeightsSum() const;

                //- Return access to normalisation factor of source
                //- patch weights (i.e. the sum before normalisation)
                inline scalarField& srcWeightsSum();

                //- Return const access to source patch face centroids
                inline const pointListList& srcCentroids() const;

                //- Return access to source patch face centroids
                inline pointListList& srcCentroids();

                //- Source map - valid only if singlePatchProc = -1
                //- This gets source data into a form to be consumed by
                //- tgtAddress, tgtWeights
                inline const mapDistribute& srcMap() const;

                //- Pointer to the source map (if distributed).
                //- Can be checked as a bool.
                inline const mapDistribute* hasSrcMap() const noexcept;


            // Target patch

                //- Return const access to target patch face areas
                inline const List<scalar>& tgtMagSf() const;

                //- Return access to target patch face areas
                inline List<scalar>& tgtMagSf();

                //- Return const access to target patch addressing
                inline const labelListList& tgtAddress() const;

                //- Return access to target patch addressing
                inline labelListList& tgtAddress();

                //- Return const access to target patch weights
                inline const scalarListList& tgtWeights() const;

                //- Return access to target patch weights
                inline scalarListList& tgtWeights();

                //- Return const access to normalisation factor of target
                //- patch weights (i.e. the sum before normalisation)
                inline const scalarField& tgtWeightsSum() const;

                //- Return access to normalisation factor of target
                //- patch weights (i.e. the sum before normalisation)
                inline scalarField& tgtWeightsSum();

                //- Target map - valid only if singlePatchProc=-1.
                //- This gets target data into a form to be consumed by
                //- srcAddress, srcWeights
                inline const mapDistribute& tgtMap() const;

                //- Pointer to the target map (if distributed).
                //- Can be checked as a bool.
                inline const mapDistribute* hasTgtMap() const noexcept;


        // Manipulation

            //- Update addressing, weights and (optional) centroids
            virtual bool calculate
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr = nullptr
            );

            //- Set the maps, addresses and weights from an external source
            void reset
            (
                autoPtr<mapDistribute>&& srcToTgtMap,
                autoPtr<mapDistribute>&& tgtToSrcMap,
                labelListList&& srcAddress,
                scalarListList&& srcWeights,
                labelListList&& tgtAddress,
                scalarListList&& tgtWeights,
                const label singlePatchProc
            );

            //- Append additional addressing and weights
            void append
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            );

            //- Normalise the weights
            void normaliseWeights(const bool conformal, const bool output);


        // Evaluation

            // Low-level

                //- Weighted sum of contributions
                template<class Type, class CombineOp>
                static void weightedSum
                (
                    const scalar lowWeightCorrection,
                    const labelListList& allSlots,
                    const scalarListList& allWeights,
                    const scalarField& weightsSum,
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues
                );

                //- Weighted sum of contributions
                template<class Type>
                void weightedSum
                (
                    const bool interpolateToSource,
                    const UList<Type>& fld,
                    List<Type>& result,
                    const UList<Type>& defaultValues
                ) const;

                //- Interpolate from target to source with supplied op
                //- to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToSource
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;

                //- Interpolate from source to target with supplied op
                //- to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToTarget
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;


            //- Interpolate from target to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target to source
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;


        // Point intersections

            //- Return source patch face index of point on target patch face
            label srcPointFace
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const vector& n,
                const label tgtFacei,
                point& tgtPoint
            )
            const;

            //- Return target patch face index of point on source patch face
            label tgtPointFace
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const vector& n,
                const label srcFacei,
                point& srcPoint
            )
            const;


        // Checks

            //- Check if src addresses are present in tgt addresses and
            //- viceversa
            bool checkSymmetricWeights(const bool log) const;

            //- Write face connectivity as OBJ file
            void writeFaceConnectivity
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const labelListList& srcAddress
            ) const;


        // I-O

            //- Write AMI as a dictionary
            virtual void write(Ostream& os) const;

            //- Write AMI raw
            virtual bool writeData(Ostream& os) const;


    // Housekeeping

        //- Non-const access to the up-to-date flag
        //  \deprecated Prefer the upToDate(bool) setter (JAN-2023)
        bool& upToDate() noexcept { return upToDate_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
