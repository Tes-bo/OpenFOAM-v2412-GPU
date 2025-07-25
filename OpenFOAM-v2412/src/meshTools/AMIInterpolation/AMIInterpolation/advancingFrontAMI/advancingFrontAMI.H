/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::advancingFrontAMI

Description
    Base class for Arbitrary Mesh Interface (AMI) methods

SourceFiles
    advancingFrontAMI.C

\*---------------------------------------------------------------------------*/

#ifndef advancingFrontAMI_H
#define advancingFrontAMI_H

#include "className.H"
#include "DynamicList.H"
#include "faceAreaIntersect.H"
#include "pointList.H"
#include "AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class advancingFrontAMI Declaration
\*---------------------------------------------------------------------------*/

class advancingFrontAMI
:
    public AMIInterpolation
{
public:

    //- Area normalisation mode
    enum class areaNormalisationMode
    {
        project,
        mag
    };

    static const Enum<areaNormalisationMode> areaNormalisationModeNames_;


private:

    // Private Member Functions

        //- No copy assignment
        void operator=(const advancingFrontAMI&) = delete;


        // Parallel operations

            label calcOverlappingProcs
            (
                const List<treeBoundBoxList>& procBb,
                const treeBoundBox& bb,
                boolList& overlaps
            ) const;

            void distributePatches
            (
                const mapDistribute& map,
                const primitivePatch& pp,
                const globalIndex& gi,
                List<faceList>& faces,
                List<pointField>& points,
                List<labelList>& tgtFaceIDs
            ) const;

            void distributeAndMergePatches
            (
                const mapDistribute& map,
                const primitivePatch& tgtPatch,
                const globalIndex& gi,
                faceList& tgtFaces,
                pointField& tgtPoints,
                labelList& tgtFaceIDs
            ) const;

            autoPtr<mapDistribute> calcProcMap
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            ) const;


protected:

    // Protected data

        //- Maximum squared distance
        const scalar maxDistance2_;

        //- Minimum (cos of) angle. 1 for perfectly matching.
        const scalar minCosAngle_;

        //- Storage for src-side triangle decomposition
        List<DynamicList<face>> srcTris_;

        //- Storage for tgt-side triangle decomposition
        List<DynamicList<face>> tgtTris_;

        //- Demand-driven extended target mesh (distributed parallel usage)
        autoPtr<primitivePatch> extendedTgtPatchPtr_;

        //- Extended patch faces
        faceList extendedTgtFaces_;

        //- Extended patch points
        pointField extendedTgtPoints_;

        //- Extended patch face IDs
        labelList extendedTgtFaceIDs_;

        //- Extended patch map
        autoPtr<mapDistribute> extendedTgtMapPtr_;

        //- Labels of faces that are not overlapped by any target faces
        //- (should be empty for correct functioning for fully covered AMIs)
        labelList srcNonOverlap_;

        //- Octree used to find face seeds
        autoPtr<indexedOctree<treeType>> treePtr_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;

        //- Area normalisation mode; default = project
        areaNormalisationMode areaNormalisationMode_;


    // Protected Member Functions

        // Helper functions

            //- Create a map that extends tgtPatch so that it covers srcPatch
            void createExtendedTgtPatch();

            //- Check AMI patch coupling
            void checkPatches() const;

            //- Is source/target a valid pair (i.e. not too far/different
            //  orientation). Used for prefiltering before e.g. area overlap
            bool isCandidate(const label srcFacei, const label tgtFacei) const;

            virtual bool calculate
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr = nullptr
            );

            //- Initialise walk and return true if all ok
            bool initialiseWalk
            (
                label& srcFacei,
                label& tgtFacei
            );

            //- Write triangle intersection to OBJ file
            void writeIntersectionOBJ
            (
                const scalar area,
                const face& f1,
                const face& f2,
                const pointField& f1Points,
                const pointField& f2Points
            ) const;


        // Common AMI method functions

            label findTargetFace
            (
                const label srcFacei,
                const UList<label>& excludeFaces = UList<label>::null(),
                const label srcFacePti = -1
            ) const;

            //- Add faces neighbouring facei to the ID list
            void appendNbrFaces
            (
                const label facei,
                const primitivePatch& patch,
                const DynamicList<label>& visitedFaces,
                DynamicList<label>& faceIDs
            ) const;

            //- Helper function to decompose a patch
            void triangulatePatch
            (
                const primitivePatch& patch,
                List<DynamicList<face>>& tris,
                List<scalar>& magSf
            ) const;

            //- Correction for non-conformal interpolations, e.g. for ACMI
            virtual void nonConformalCorrection();


public:

    //- Runtime type information
    TypeName("advancingFrontAMI");

    // Constructors

        //- Construct from components
        advancingFrontAMI
        (
            const dictionary& dict,
            const bool reverseTarget
        );

        //- Construct from components
        advancingFrontAMI
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh
        );

        //- Construct as copy
        advancingFrontAMI(const advancingFrontAMI& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>(new advancingFrontAMI(*this));
        }


    //- Destructor
    virtual ~advancingFrontAMI() = default;


    // Member Functions

        //- Return const access to the source patch
        inline const primitivePatch& srcPatch() const;

        //- Return const access to the target patch
        inline const primitivePatch& tgtPatch() const;

        //- Labels of faces that are not overlapped by any target faces
        //  Note: this should be empty for correct functioning
        inline const labelList& srcNonOverlap() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "advancingFrontAMII.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
