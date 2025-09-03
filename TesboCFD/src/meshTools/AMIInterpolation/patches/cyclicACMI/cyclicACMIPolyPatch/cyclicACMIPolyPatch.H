/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::cyclicACMIPolyPatch

Description
    Cyclic patch for Arbitrarily Coupled Mesh Interface (ACMI).

    Mixes cyclicAMI behaviour with non-coupled patch behaviour using
    the overlap area fraction. The non-coupled patch is specified through
    the nonOverlapPatch keyword.

Usage
    Example of the patch specification:

        type            cyclicACMI;
        neighbourPatch  ACMI2_couple;       // cyclicAMI neighbour patch
        nonOverlapPatch ACMI1_blockage;     // patch for uncoupled faces

        // Optional time-dependent scaling (PatchFunction1)
        scale           table
        (
            (0.00   1.0)
            (0.02   1.0)
            (0.0201 0.0)
        );

See also
    cyclicAMIPolyPatch.C

SourceFiles
    cyclicACMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIPolyPatch_H
#define cyclicACMIPolyPatch_H

#include "cyclicAMIPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"
#include "PatchFunction1.H"
#include "uniformDimensionedFields.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cyclicACMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIPolyPatch
:
    public cyclicAMIPolyPatch
{

private:

    // Private data

        //- Name of non-overlapping patch
        const word nonOverlapPatchName_;

        //- Index of non-overlapping patch
        mutable label nonOverlapPatchID_;

        //- Mask/weighting for source patch
        mutable scalarField srcMask_;

        //- Mask/weighting for target patch
        mutable scalarField tgtMask_;


        // Scaling of overlap (optional)

            //- Weighting for source mask
            mutable autoPtr<PatchFunction1<scalar>> srcScalePtr_;

            //- Weighting for target mask
            mutable autoPtr<PatchFunction1<scalar>> tgtScalePtr_;

            //- Stored raw/non-scaled face areas
            mutable vectorField thisSf_;
            mutable vectorField thisNoSf_;
            mutable vectorField nbrSf_;
            mutable vectorField nbrNoSf_;

            //- Scaled version of source mask
            mutable scalarField srcScaledMask_;

            //- Scaled version of target mask
            mutable scalarField tgtScaledMask_;

            //- Flag to detect whether AMI is up to date with mesh points
            mutable uniformDimensionedScalarField AMITime_;

            //- Flag to detect whether scaled masks are up to date with
            //  current time
            mutable label prevTimeIndex_;


protected:

    // Protected Member Functions

        //- Report coverage statics, e.g. number of uncovered/blended/covered
        //- faces
        void reportCoverage
        (
            const word& name,
            const scalarField& weightSum
        ) const;

        //- Reset the AMI interpolator, supply patch points
        virtual void resetAMI(const UList<point>& points) const;

        //-  Reset the AMI interpolator, use current patch points
        virtual void resetAMI() const;

        //- Scale patch face areas to maintain physical area
        virtual void scalePatchFaceAreas();

        //- Scale patch face areas to maintain physical area
        virtual void scalePatchFaceAreas
        (
            const cyclicACMIPolyPatch& acmipp,
            const scalarField& mask,
            const vectorList& faceArea,
            const vectorList& noFaceArea
        );

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

        //- Clear geometry
        virtual void clearGeom();


public:

    //- Runtime type information
    TypeName("cyclicACMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicACMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct from dictionary
        cyclicACMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName,
            const word& nonOverlapPatchName
        );

        //- Construct given the original patch and a map
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<cyclicACMIPolyPatch>(*this, bm);
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicACMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    neighbPatchName(),
                    nonOverlapPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicACMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicACMIPolyPatch() = default;


    // Member Functions

        // Implicit treatment functions

            //- Return number of new internal sub-faces and new proc faces
            virtual void newInternalProcFaces(label&, label&) const;

            //- Return collocated faces
            virtual refPtr<labelListList> mapCollocatedFaces() const;


        // Access

            //- Return a reference to the neighbour patch
            virtual const cyclicACMIPolyPatch& neighbPatch() const;

            //- Non-overlapping patch name
            inline const word& nonOverlapPatchName() const;

            //- Non-overlapping patch ID
            virtual label nonOverlapPatchID() const;

            //- Return a const reference to the non-overlapping patch
            inline const polyPatch& nonOverlapPatch() const;

            //- Return a reference to the non-overlapping patch
            inline polyPatch& nonOverlapPatch();

            //- Mask field where 1 = overlap(coupled), 0 = no-overlap
            inline const scalarField& mask() const;

            //- Return the mask/weighting for the source patch
            virtual const scalarField& srcMask() const;

            //- Return the mask/weighting for the target patch
            virtual const scalarField& tgtMask() const;

            //- Overlap tolerance
            inline static scalar tolerance();


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;

        // Handling optional scaling (time dependency)

            //- Update the AMI and patch areas. Return true if anything
            //  updated
            virtual bool updateAreas() const;

            //- Return true if given object is up to date with *this
            //  (note : like regIOobject::upToDate but operates on object)
            bool upToDate(const regIOobject&) const;

            //- Set object up to date with *this
            //  (note : like regIOobject::setUpToDate but operates on object)
            void setUpToDate(regIOobject&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cyclicACMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
