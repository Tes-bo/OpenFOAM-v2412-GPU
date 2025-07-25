/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::particle

Description
    Base particle class

\*---------------------------------------------------------------------------*/

#ifndef Foam_particle_H
#define Foam_particle_H

#include "vector.H"
#include "barycentric.H"
#include "barycentricTensor.H"
#include "Cloud.H"
#include "IDLList.H"
#include "pointField.H"
#include "faceList.H"
#include "OFstream.H"
#include "FixedList.H"
#include "polyMeshTetDecomposition.H"
#include "particleMacros.H"
#include "vectorTensorTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class particle;
class polyPatch;
class cyclicPolyPatch;
class cyclicAMIPolyPatch;
class cyclicACMIPolyPatch;
class processorPolyPatch;
class symmetryPlanePolyPatch;
class symmetryPolyPatch;
class wallPolyPatch;
class wedgePolyPatch;

Ostream& operator<<(Ostream&, const particle&);
bool operator==(const particle&, const particle&);
bool operator!=(const particle&, const particle&);

/*---------------------------------------------------------------------------*\
                          Class Particle Declaration
\*---------------------------------------------------------------------------*/

class particle
:
    public IDLList<particle>::link
{
    // Private Data

        //- Size in bytes of the position data
        static const std::size_t sizeofPosition;

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;

        //- The value of nBehind_ at which tracking is abandoned. See the
        //  description of nBehind_.
        static const label maxNBehind_;


public:

    class trackingData
    {
    public:

        // Public Data

            //- Flag to switch processor
            bool switchProcessor;

            //- Flag to indicate whether to keep particle (false = delete)
            bool keepParticle;


        // Constructor
        template<class TrackCloudType>
        trackingData(const TrackCloudType& cloud)
        {}
    };


    //- Old particle positions content for OpenFOAM-1706 and earlier
    struct positionsCompat1706
    {
        vector position;
        label celli;
        label facei;
        scalar stepFraction;
        label tetFacei;
        label tetPti;
        label origProc;
        label origId;
    };


private:

    // Private Data

        //- Reference to the polyMesh database
        const polyMesh& mesh_;

        //- Coordinates of particle
        barycentric coordinates_;

        //- Index of the cell it is in
        label celli_;

        //- Index of the face that owns the decomposed tet that the
        //- particle is in
        label tetFacei_;

        //- Index of the point on the face that defines the decomposed
        //- tet that the particle is in.  Relative to the face base
        //- point.
        label tetPti_;

        //- Face index if the particle is on a face otherwise -1
        label facei_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- The distance behind the maximum distance reached so far
        scalar behind_;

        //- The number of tracks carried out that ended in a distance behind the
        //- maximum distance reached so far. Once this reaches maxNBehind_,
        //  tracking is abandoned for the current step. This is needed because
        //  when tetrahedra are inverted a straight trajectory can form a closed
        //  loop through regions of overlapping positive and negative space.
        //  Without this break clause, such loops can result in a valid track
        //  which never ends.
        label nBehind_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


    // Private Member Functions

        // Tetrahedra functions

            //- Get the transformation associated with the current tet. This
            //  will convert a barycentric position within the tet to a
            //  Cartesian position in the global coordinate system. The
            //  conversion is x = A & y, where x is the Cartesian position, y is
            //  the barycentric position and A is the transformation tensor.
            inline barycentricTensor stationaryTetTransform() const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = (x - centre) & T. The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the transposed inverse of the forward transform tensor, A,
            //  multiplied by its determinant, detA. This separation allows
            //  the barycentric tracking algorithm to function on inverted or
            //  degenerate tetrahedra.
            void stationaryTetReverseTransform
            (
                vector& centre,
                scalar& detA,
                barycentricTensor& T
            ) const;

            //- Get the vertices of the current moving tet. Two values are
            //  returned for each vertex. The first is a constant, and the
            //  second is a linear coefficient of the track fraction.
            inline void movingTetGeometry
            (
                const scalar endStepFraction,
                Pair<vector>& centre,
                Pair<vector>& base,
                Pair<vector>& vertex1,
                Pair<vector>& vertex2
            ) const;

            //- Get the transformation associated with the current, moving, tet.
            //  This is of the same form as for the static case. As with the
            //  moving geometry, a linear function of the tracking fraction is
            //  returned for each component.
            inline Pair<barycentricTensor> movingTetTransform
            (
                const scalar endStepFraction
            ) const;

            //- Get the reverse transformation associated with the current,
            //  moving, tet. This is of the same form as for the static case. As
            //  with the moving geometry, a function of the tracking fraction is
            //  returned for each component. The functions are higher order than
            //  for the forward transform; the determinant is cubic, and the
            //  tensor is quadratic.
            void movingTetReverseTransform
            (
                const scalar endStepFraction,
                Pair<vector>& centre,
                FixedList<scalar, 4>& detA,
                FixedList<barycentricTensor, 3>& T
            ) const;


        // Transformations

            //- Reflection transform. Corrects the coordinates when the particle
            //  moves between two tets which share a base vertex, but for which
            //  the other two non cell-centre vertices are reversed. All hits
            //  which retain the same face behave this way, as do face hits.
            void reflect();

            //- Rotation transform. Corrects the coordinates when the particle
            //  moves between two tets with different base vertices, but are
            //  otherwise similarly oriented. Hits which change the face within
            //  the cell make use of both this and the reflect transform.
            void rotate(const bool direction);


        // Topology changes

            //- Change tet within a cell. Called after a triangle is hit.
            void changeTet(const label tetTriI);

            //- Change tet face within a cell. Called by changeTet.
            void changeFace(const label tetTriI);

            //- Change cell. Called when the particle hits an internal face.
            void changeCell();

            //- Put the particle on the lowest indexed patch for the current set
            //  of coincident faces. In the case of an ACMI-wall pair, this will
            //  move the particle from the wall face to the ACMI face, because
            //  ACMI patches are always listed before their associated non-
            //  overlapping patch.
            void changeToMasterPatch();


        // Geometry changes

            //- Locate the particle at the given position
            void locate
            (
                const vector& position,
                const vector* direction,
                const label celli,
                const bool boundaryFail,
                const string& boundaryMsg
            );


protected:

    // Patch Interactions

        //- Read particle from stream. Optionally (for old format) return
        //  read position. Used by construct-from-Istream
        void readData
        (
            Istream& is,
            point& position,
            const bool readFields,
            const bool newFormat,
            const bool doLocate
        );

        //- Overridable function to handle the particle hitting a patch.
        //  Executed before other patch-hitting functions.
        template<class TrackCloudType>
        bool hitPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a wedgePatch
        template<class TrackCloudType>
        void hitWedgePatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a
        //  symmetryPlanePatch
        template<class TrackCloudType>
        void hitSymmetryPlanePatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a symmetryPatch
        template<class TrackCloudType>
        void hitSymmetryPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackCloudType>
        void hitCyclicPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a cyclicAMIPatch
        template<class TrackCloudType>
        void hitCyclicAMIPatch(TrackCloudType&, trackingData&, const vector&);

        //- Overridable function to handle the particle hitting a
        //  cyclicACMIPatch
        template<class TrackCloudType>
        void hitCyclicACMIPatch(TrackCloudType&, trackingData&, const vector&);

        //- Overridable function to handle the particle hitting a processorPatch
        template<class TrackCloudType>
        void hitProcessorPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackCloudType>
        void hitWallPatch(TrackCloudType&, trackingData&);


        //- Dispatch function for boundary face interaction. Calls one of
        //  the above (hitWedgePatch, hitCyclicPatch etc) depending on the
        //  patch type
        template<class TrackCloudType>
        void hitBoundaryFace
        (
            const vector& direction,
            TrackCloudType& cloud,
            trackingData& td
        );


public:

    // Static Data Members

        //- Runtime type information
        TypeName("particle");

        //- String representation of properties
        DefinePropertyList
        (
            "(coordinatesa coordinatesb coordinatesc coordinatesd) "
            "celli tetFacei tetPti facei stepFraction origProc origId"
        );

        //- Cumulative particle counter - used to provide unique ID
        static label particleCount_;

        //- Write particle coordinates file (v1712 and later)
        //- Default is true
        static bool writeLagrangianCoordinates;

        //- Write particle positions file (v1706 format and earlier)
        //- Default is true (disable in etc/controlDict)
        static bool writeLagrangianPositions;


    // Constructors

        //- Construct from components
        particle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell.
        //  Searches for the rest of the required topology.
        particle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        );

        //- Construct from position components
        particle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const bool doLocate = true
        );


        //- Construct from Istream
        particle
        (
            const polyMesh& mesh,
            Istream&,
            const bool readFields = true,
            const bool newFormat = true,
            const bool doLocate = true
        );

        //- Construct as a copy with reference to a mesh
        particle(const particle& p, const polyMesh& mesh);

        //- Copy construct
        particle(const particle& p);

        //- Construct a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>::New(*this);
        }


    // Factory Methods

        //- Clone a particle
        template<class Derived>
        static autoPtr<particle> Clone(const Derived& p)
        {
            return autoPtr<particle>(new Derived(p));
        }

        //- Clone a particle with a mesh reference
        template<class Derived>
        static autoPtr<particle> Clone(const Derived& p, const polyMesh& mesh)
        {
            return autoPtr<particle>(new Derived(p, mesh));
        }

        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<particle> operator()(Istream& is) const
            {
                return autoPtr<particle>::New(mesh_, is, true);
            }
        };


    //- Destructor
    virtual ~particle() = default;


    // Member Functions

        // Access

            //- Get unique particle creation id
            inline label getNewParticleID() const;

            //- Return the mesh database
            inline const polyMesh& mesh() const noexcept;

            //- Return current particle coordinates
            inline const barycentric& coordinates() const noexcept;

            //- Return current cell particle is in
            inline label cell() const noexcept;

            //- Return current cell particle is in for manipulation
            inline label& cell() noexcept;

            //- Return current tet face particle is in
            inline label tetFace() const noexcept;

            //- Return current tet face particle is in for manipulation
            inline label& tetFace() noexcept;

            //- Return current tet face particle is in
            inline label tetPt() const noexcept;

            //- Return current tet face particle is in for manipulation
            inline label& tetPt() noexcept;

            //- Return current face particle is on otherwise -1
            inline label face() const noexcept;

            //- Return current face particle is on for manipulation
            inline label& face() noexcept;

            //- Return the fraction of time-step completed
            inline scalar stepFraction() const noexcept;

            //- Return the fraction of time-step completed
            inline scalar& stepFraction() noexcept;

            //- Return the originating processor ID
            inline label origProc() const noexcept;

            //- Return the originating processor ID
            inline label& origProc() noexcept;

            //- Return the particle ID on the originating processor
            inline label origId() const noexcept;

            //- Return the particle ID on the originating processor
            inline label& origId() noexcept;


        // Check

            //- Return the step fraction change within the overall time-step.
            //  Returns the start value and the change as a scalar pair. Always
            //  return Pair<scalar>(0, 1), unless sub-cycling is in effect, in
            //  which case the values will reflect the span of the sub-cycle
            //  within the time-step.
            inline Pair<scalar> stepFractionSpan() const;

            //- Return the current fraction within the timestep. This differs
            //  from the stored step fraction due to sub-cycling.
            inline scalar currentTimeFraction() const;

            //- Return indices of the current tet that the particle occupies.
            inline tetIndices currentTetIndices() const noexcept;

            //- Return the current tet transformation tensor
            inline barycentricTensor currentTetTransform() const;

            //- The (unit) normal of the tri on tetFacei_ for the current tet.
            inline vector normal() const;

            //- Is the particle on a face?
            inline bool onFace() const noexcept;

            //- Is the particle on an internal face?
            inline bool onInternalFace() const noexcept;

            //- Is the particle on a boundary face?
            inline bool onBoundaryFace() const noexcept;

            //- Return the index of patch that the particle is on
            inline label patch() const;

            //- Return current particle position
            inline vector position() const;

            //- Reset particle data
            inline void reset();


    // Track

        //- Track along the displacement for a given fraction of the overall
        //  step. End when the track is complete, or when a boundary is hit.
        //  On exit, stepFraction_ will have been incremented to the current
        //  position, and facei_ will be set to the index of the boundary
        //  face that was hit, or -1 if the track completed within a cell.
        //  The proportion of the displacement still to be completed is
        //  returned.
        scalar track
        (
            const vector& displacement,
            const scalar fraction
        );

        //- As particle::track, but also stops on internal faces.
        scalar trackToFace
        (
            const vector& displacement,
            const scalar fraction
        );

        //- As particle::trackToFace, but also stops on tet triangles. On
        //  exit, tetTriI is set to the index of the tet triangle that was
        //  hit, or -1 if the end position was reached.
        scalar trackToTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- As particle::trackToTri, but for stationary meshes
        scalar trackToStationaryTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- As particle::trackToTri, but for moving meshes
        scalar trackToMovingTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- Hit the current face. If the current face is internal than this
        //  crosses into the next cell. If it is a boundary face then this will
        //  interact the particle with the relevant patch.
        template<class TrackCloudType>
        void hitFace
        (
            const vector& direction,
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Convenience function. Combines trackToFace and hitFace
        template<class TrackCloudType>
        void trackToAndHitFace
        (
            const vector& direction,
            const scalar fraction,
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Get the displacement from the mesh centre. Used to correct the
        //  particle position in cases with reduced dimensionality. Returns a
        //  zero vector for three-dimensional cases.
        vector deviationFromMeshCentre() const;


    // Patch data

        //- Get the normal and velocity of the current patch location
        inline void patchData(vector& n, vector& U) const;


    // Transformations

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);


    // Parallel transfer

        //- Convert global addressing to the processor patch local equivalents
        void prepareForParallelTransfer();

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        void correctAfterParallelTransfer(const label patchi, trackingData& td);


    // Interaction list referral

        //- Break the topology and store the particle position so that the
        //  particle can be referred.
        void prepareForInteractionListReferral
        (
            const vectorTensorTransform& transform
        );

        //- Correct the topology after referral. The particle may still be
        //  outside the stored tet and therefore not track-able.
        void correctAfterInteractionListReferral(const label celli);


    // Decompose and reconstruct

        //- Return the tet point appropriate for decomposition or reconstruction
        //  to or from the given mesh.
        label procTetPt
        (
            const polyMesh& procMesh,
            const label procCell,
            const label procTetFace
        ) const;


    // Mapping

        //- Map after a topology change
        void autoMap(const vector& position, const mapPolyMesh& mapper);

        //- Set the addressing based on the provided position
        void relocate(const point& position, const label celli = -1);



    // I-O

        //- Write the name representation to stream
        template<class Type>
        static void writePropertyName
        (
            Ostream& os,
            const word& name,
            const word& delim
        );

        //- Write a named particle property to stream,
        //- optionally filtered based on its name
        template<class Type>
        static void writeProperty
        (
            Ostream& os,
            const word& name,
            const Type& value,
            const bool nameOnly,
            const word& delim,
            const wordRes& filters = wordRes::null()
        );

        //- Write a named particle property list to stream,
        //- optionally filtered based on its name
        template<class Type>
        static void writeProperty
        (
            Ostream& os,
            const word& name,
            const Field<Type>& values,
            const bool nameOnly,
            const word& delim,
            const wordRes& filters = wordRes::null()
        );

        //- Read the fields associated with the owner cloud
        template<class TrackCloudType>
        static void readFields(TrackCloudType& c);

        //- Write the fields associated with the owner cloud
        template<class TrackCloudType>
        static void writeFields(const TrackCloudType& c);

        //- Write individual particle properties to stream
        void writeProperties
        (
            Ostream& os,
            const wordRes& filters,
            const word& delim,
            const bool namesOnly
        ) const;

        //- Read particle fields as objects from the obr registry
        template<class CloudType>
        static void readObjects(CloudType& c, const objectRegistry& obr);

        //- Write particle fields as objects into the obr registry
        //  Always writes "position", not "coordinate"
        template<class CloudType>
        static void writeObjects(const CloudType& c, objectRegistry& obr);

        //- Write the particle barycentric coordinates and cell info
        void writeCoordinates(Ostream& os) const;

        //- Write the particle position and cell id
        virtual void writePosition(Ostream& os) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const particle&);

        friend bool operator==(const particle& pA, const particle& pB);

        friend bool operator!=(const particle& pA, const particle& pB);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
