/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::cyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

    Includes provision for updating the patch topology to enforce a 1-to-1
    face match across the interface, based on the \c createAMIFaces flag.

    The manipulations are based on the reference:

    \verbatim
    H.J. Aguerre, S. Márquez Damián, J.M. Gimenez, N.M.Nigro, Conservative
    handling of arbitrary non-conformal interfaces using an efficient
    supermesh, Journal of Computational Physics 335(15) 21-49. 2017.
    https://doi.org/10.1016/j.jcp.2017.01.018.
    \endverbatim

SourceFiles
    cyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicAMIPolyPatch_H
#define Foam_cyclicAMIPolyPatch_H

#include "AMIInterpolation.H"
#include "coupledPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"
#include "coupleGroupIdentifier.H"
#include "faceAreaWeightAMI.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPolyPatch
:
    public coupledPolyPatch
{
    // Private Member Functions

        //- Return normal of face at max distance from rotation axis
        vector findFaceNormalMaxRadius(const pointField& faceCentres) const;

        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );


protected:

    // Protected Data

        //- Name of other half
        mutable word nbrPatchName_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable label nbrPatchID_;

        //- Particle displacement fraction accross AMI
        const scalar fraction_;


        // Transformations

            // For rotation

                //- Axis of rotation for rotational cyclics
                vector rotationAxis_;

                //- Point on axis of rotation for rotational cyclics
                point rotationCentre_;

                //- Flag to show whether the rotation angle is defined
                bool rotationAngleDefined_;

                //- Rotation angle
                scalar rotationAngle_;


            // For translation

                //- Translation vector
                vector separationVector_;


            // Triggering periodic interpolation

                //- Periodic patch name
                word periodicPatchName_;

                //- Periodic patch
                mutable label periodicPatchID_;


        //- AMI interpolation class
        mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;

        //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;


        // Change of topology as AMI is updated

            //- Flag to indicate that new AMI faces will created
            //  Set by the call to changeTopology
            mutable bool createAMIFaces_;

            //- Move face centres (default = no)
            bool moveFaceCentres_;

            mutable bool updatingAMI_;

            labelListList srcFaceIDs_;

            labelListList tgtFaceIDs_;

            //- Temporary storage for AMI face areas
            mutable vectorField faceAreas0_;

            //- Temporary storage for AMI face centres
            mutable vectorField faceCentres0_;


    // Protected Member Functions

        // Topology change

            //- Collect faces to remove in the topoChange container
            virtual bool removeAMIFaces(polyTopoChange& topoChange);

            //- Collect faces to add in the topoChange container
            virtual bool addAMIFaces(polyTopoChange& topoChange);

            //- Set properties of newly inserted faces after topological changes
            virtual void setAMIFaces();

            //- Helper to re-apply the geometric scaling lost during mesh
            //- updates
            virtual void restoreScaledGeometry();


        //- Create and return pointer to the projection surface
        const autoPtr<searchableSurface>& surfPtr() const;

        //- Create a coordinate system from the periodic patch (or nullptr)
        autoPtr<coordSystem::cylindrical> cylindricalCS() const;

        //- Reset the AMI interpolator, supply patch points
        virtual void resetAMI(const UList<point>& points) const;

        //-  Reset the AMI interpolator, use current patch points
        virtual void resetAMI() const;

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

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
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicAMIPolyPatch
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
        cyclicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicAMIPolyPatch(const cyclicAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<cyclicAMIPolyPatch>(*this, bm);
        }

        //- Construct and return a clone, resetting the face list
        //- and boundary mesh
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
                new cyclicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
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
                new cyclicAMIPolyPatch
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
    virtual ~cyclicAMIPolyPatch() = default;


    // Member Functions

        // Implicit treatment functions

            //- Return number of new internal of this polyPatch faces
            virtual void newInternalProcFaces(label&, label&) const;


            //- Return nbrCells
            virtual const labelUList& nbrCells() const
            {
                return neighbPatch().faceCells();
            }

            //- Return nbr patch ID
            virtual label neighbPolyPatchID() const
            {
                return neighbPatch().index();
            }

            //- Return collocated faces map
            virtual refPtr<labelListList> mapCollocatedFaces() const
            {
                const labelListList& sourceFaces = AMI().srcAddress();
                return refPtr<labelListList>(new labelListList(sourceFaces));
            }

            //- Return implicit master
            virtual bool masterImplicit() const
            {
                return owner();
            }


        // Access

            //- Tolerance used e.g. for area calculations/limits
            static const scalar tolerance_;

            //- Flag to indicate whether the AMI can be reset
            inline bool canResetAMI() const;

            //- Return access to the createAMIFaces flag
            inline bool createAMIFaces() const;

            //- Return access to the updated flag
            inline bool updatingAMI() const;

            //- Return true if this patch changes the mesh topology
            //  True when createAMIFaces is true
            virtual bool changeTopology() const;

            //- Set topology changes in the polyTopoChange object
            virtual bool setTopology(polyTopoChange& topoChange);

            //- Is patch 'coupled'. Note that on AMI the geometry is not
            //- coupled but the fields are!
            virtual bool coupled() const
            {
                return false;
            }

            //- Neighbour patch name
            inline const word& neighbPatchName() const;

            //- Neighbour patch ID
            virtual label neighbPatchID() const;

            //- Particle fraction increase between AMI pathces
            inline scalar fraction() const;

            //- Does this side own the patch?
            virtual bool owner() const;

            //- Return a reference to the neighbour patch
            virtual const cyclicAMIPolyPatch& neighbPatch() const;

            //- Periodic patch ID (or -1)
            label periodicPatchID() const;

            //- Return a reference to the AMI interpolator
            const AMIPatchToPatchInterpolation& AMI() const;

            //- Helper function to return the weights
            inline const scalarListList& weights() const;

            //- Helper function to return the weights sum
            inline const scalarField& weightsSum() const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;

            //- Return access to the initial face areas
            //  Used for topology change
            inline vectorField& faceAreas0() const;

            //- Return access to the initial face centres
            //  Used for topology change
            inline vectorField& faceCentres0() const;


            // Transformations

                //- Axis of rotation for rotational cyclic AMI
                inline const vector& rotationAxis() const;

                //- Point on axis of rotation for rotational cyclic AMI
                inline const point& rotationCentre() const;

                //- Translation vector for translational cyclic AMI
                inline const vector& separationVector() const;

                //- Transform patch-based positions from nbr side to this side
                virtual void transformPosition(pointField&) const;

                //- Transform a patch-based position from nbr side to this side
                virtual void transformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based position from this side to nbr side
                virtual void reverseTransformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based direction from this side to
                //- nbr side
                virtual void reverseTransformDirection
                (
                    vector& d,
                    const label facei
                ) const;


            // Interpolations

                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate without periodic
                template<class Type>
                tmp<Field<Type>> interpolateUntransformed
                (
                    const Field<Type>& fld,
                    const UList<Type>& defaultValues
                ) const;

                //- Low-level interpolate List
                template<class Type, class CombineOp>
                void interpolate
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;


            // Interpolations (non-blocking). Subject to change

                template<class Type>
                void initInterpolateUntransformed
                (
                    const Field<Type>& fld,
                    labelRange& sendRequests,
                    PtrList<List<Type>>& sendBuffers,
                    labelRange& recvRequests,
                    PtrList<List<Type>>& recvBuffers
                ) const;

                template<class Type>
                void initInterpolate
                (
                    const Field<Type>& fld,
                    labelRange& sendRequests,
                    PtrList<List<Type>>& sendBuffers,
                    labelRange& recvRequests,
                    PtrList<List<Type>>& recvBuffers
                ) const;

                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& localFld,
                    const labelRange& requests,  // The receive requests
                    const PtrList<List<Type>>& recvBuffers,
                    const UList<Type>& defaultValues
                ) const;


        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialize ordering for primitivePatch. Does not
        //- refer to *this (except for name() and type() etc.)
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

        //- Return face index on neighbour patch which shares point p
        //- following trajectory vector n
        label pointFace
        (
            const label facei,
            const vector& n,
            point& p
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cyclicAMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicAMIPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
