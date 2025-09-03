/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::mappedPatchBase

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

    If constructed from dictionary:
    \verbatim
        // Optional world to sample (default is all)
        //sampleWorld solidSim;

        // Optional explicit coupling (requires functionObject to synchronise
        // databases. Default is close coupling (bc to bc)
        //sampleDatabase true;

        // Region to sample (default is region0)
        sampleRegion region0;

        // What to sample:
        // - nearestCell         : sample cell containing point
        // - nearestOnlyCell     : nearest sample cell (even if not containing
        //                         point)
        // - nearestPatchFace    : nearest face on selected patch
        // - nearestPatchFaceAMI : nearest face on selected patch
                                   - patches need not conform
                                   - uses AMI interpolation
        // - nearestFace         : nearest boundary face on any patch
        // - nearestPatchPoint   : nearest patch point (for coupled points
        //                         this might be any of the points so you have
        //                         to guarantee the point data is synchronised
        //                         beforehand)
        sampleMode nearestCell;

        // If sampleMode is nearestPatchFace : patch to find faces of
        samplePatch movingWall;

        // If sampleMode is nearestPatchFace : specify patchgroup to find
        // samplePatch and sampleRegion (if not provided)
        coupleGroup baffleGroup;

        // How to supply offset (w.r.t. my patch face centres):
        // - uniform : single offset vector
        // - nonuniform : per-face offset vector
        // - normal : using supplied distance and face normal
        offsetMode uniform;

        // According to offsetMode (see above) supply one of
        // offset, offsets or distance
        offset  (1 0 0);
    \endverbatim

    Note: if offsetMode is \c normal it uses outwards pointing normals. So
    supply a negative distance if sampling inside the domain.

Note
    Storage is not optimal. It temporary collects all (patch)face centres
    on all processors to keep the addressing calculation simple.

SourceFiles
    mappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBase_H
#define mappedPatchBase_H

#include "pointField.H"
#include "Tuple2.H"
#include "pointIndexHit.H"
#include "AMIPatchToPatchInterpolation.H"
#include "coupleGroupIdentifier.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchBase
{
public:

    // Type enumerations

        //- Mesh items to sample
        enum sampleMode
        {
            NEARESTCELL,         //!< nearest cell containing sample
            NEARESTPATCHFACE,    //!< nearest face on selected patch
            NEARESTPATCHFACEAMI, //!< nearest patch face + AMI interpolation
            NEARESTPATCHPOINT,   //!< nearest point on selected patch
            NEARESTFACE,         //!< nearest face
            NEARESTONLYCELL      //!< nearest cell (even if not containing cell)
        };

        //- How to project face centres
        enum offsetMode
        {
            UNIFORM,            //!< single offset vector
            NONUNIFORM,         //!< per-face offset vector
            NORMAL              //!< use face normal + distance
        };

        static const Enum<sampleMode> sampleModeNames_;

        static const Enum<offsetMode> offsetModeNames_;


    //- Helper class for finding nearest
    //  Nearest:
    //  - point+local index
    //  - sqr(distance)
    //  - processor
    typedef Tuple2<pointIndexHit, Tuple2<scalar, label>> nearInfo;

    class nearestEqOp
    {
    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second().first() < x.second().first())
                {
                    x = y;
                }
            }
        }
    };

    class maxProcEqOp
    {
    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second().second() > x.second().second())
                {
                    x = y;
                }
            }
        }
    };


    //- nearest + world
    //  Used to only look at entries from same world
    typedef Tuple2<nearInfo, label> nearInfoWorld;

    class nearestWorldEqOp
    {
    public:

        void operator()(nearInfoWorld& x, const nearInfoWorld& y) const
        {
            // Is there a hit and is it sampling the same world
            const nearInfo& xi = x.first();
            const nearInfo& yi = y.first();
            if (yi.first().hit())
            {
                if (x.second() == y.second())
                {
                    if (!xi.first().hit())
                    {
                        x = y;
                    }
                    else if (yi.second().first() < xi.second().first())
                    {
                        x = y;
                    }
                }
            }
        }
    };


protected:

    // Protected Data

        //- Patch to sample
        const polyPatch& patch_;

        //- World to sample
        mutable word sampleWorld_;

        //- Region to sample
        mutable word sampleRegion_;

        //- What to sample
        const sampleMode mode_;

        //- Patch (if in sampleMode NEARESTPATCH*)
        mutable word samplePatch_;

        //- PatchGroup (if in sampleMode NEARESTPATCH*)
        const coupleGroupIdentifier coupleGroup_;

        //- Empty or location of database
        const autoPtr<fileName> sampleDatabasePtr_;

        //- How to obtain samples
        offsetMode offsetMode_;

        //- Offset vector (uniform)
        vector offset_;

        //- Offset vector (nonuniform)
        vectorField offsets_;

        //- Offset distance (normal)
        scalar distance_;

        //- Communicator
        label communicator_;

        //- Same region
        mutable bool sameRegion_;


        // Derived information

            //- Communication schedule:
            //
            //    - Cells/faces to sample per processor
            //    - Patch faces to receive per processor
            //    - schedule
            mutable autoPtr<mapDistribute> mapPtr_;


        // AMI interpolator (only for NEARESTPATCHFACEAMI)

            //- Flag to indicate that slave patch should be reversed for AMI
            const bool AMIReverse_;

            //- Pointer to AMI interpolator
            mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;

            //- Pointer to projection surface employed by AMI interpolator
            mutable autoPtr<searchableSurface> surfPtr_;

            //- Dictionary storing projection surface description
            dictionary surfDict_;


        // Mesh update IOobjects

            //- Local mesh update time
            mutable autoPtr<uniformDimensionedScalarField> updateMeshTimePtr_;

            //- Sample mesh update time
            mutable autoPtr<uniformDimensionedScalarField>
                updateSampleMeshTimePtr_;


    // Protected Member Functions

        // polyPatch callbacks

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry(PstreamBuffers&)
            {}

            //- Calculate the patch geometry
            virtual void calcGeometry(PstreamBuffers&)
            {}

            //- Initialise the patches for moving points
            virtual void initMovePoints(PstreamBuffers&, const pointField&)
            {}

            //- Correct patches after moving points
            virtual void movePoints(PstreamBuffers&, const pointField&)
            {}

            //- Initialise the update of the patch topology
            virtual void initUpdateMesh(PstreamBuffers&)
            {}

            //- Update of the patch topology
            virtual void updateMesh(PstreamBuffers&);


        //- Add a world-world connection
        bool addWorldConnection();

        //- Get the communicator for the world-world connection
        label getWorldCommunicator() const;

        //- Lookup mesh
        const polyMesh& lookupMesh(const word& region) const;

        //- Lookup patch
        const polyPatch& lookupPatch
        (
            const word& sampleRegion,
            const word& samplePatch
        ) const;


        //- Get the points from face-centre-decomposition face centres
        //- and project them onto the face-diagonal-decomposition triangles.
        tmp<pointField> facePoints(const polyPatch&) const;

        //- Collect single list of samples and originating processor+face +
        //- wanted world
        void collectSamples
        (
            const label mySampleWorld,  // My wanted sampling world
            const pointField& facePoints,
            pointField& samples,        // Per sample: coordinate
            labelList& patchFaceWorlds, // Per sample: wanted world
            labelList& patchFaceProcs,  // Per sample: originating proc
            labelList& patchFaces,      // Per sample: originating face
            pointField& patchFc
        ) const;

        //- Find (local) cells/faces containing samples
        void findLocalSamples
        (
            const sampleMode mode,
            const label sampleWorld,        // my world as index
            const word& sampleRegion,
            const word& samplePatch,
            const pointField& samplePoints,

            List<nearInfoWorld>& nearest
        ) const;

        //- Find (global) cells/faces containing samples
        void findSamples
        (
            const sampleMode mode,      // search mode
            const label myWorldIndex,   // my world (in index form)

            const pointField&,
            const labelList& wantedWorlds,
            const labelList& origProcs, // per sample the originating proc

            labelList& sampleProcs,     // processor containing sample
            labelList& sampleIndices,   // local index of cell/face
            pointField& sampleLocations // actual representative location
        ) const;

        //- Get the sample points given the face points
        tmp<pointField> samplePoints(const pointField&) const;

        //- Calculate mapping
        void calcMapping() const;

        //- Calculate AMI interpolator
        void calcAMI() const;


        // Database Handling

            //- Read optional database name from dictionary
            static autoPtr<fileName> readDatabase(const dictionary& dict);

            //- Lookup (sub)objectRegistry by following names of sub registries.
            //- Creates non-existing intermediate ones.
            static const objectRegistry& subRegistry
            (
                const objectRegistry& obr,
                const wordList& names,
                const label index
            );

            //- Attempt to detect an IOField<Type> and write to dictionary
            template<class Type>
            static bool writeIOField
            (
                const regIOobject& obj,
                dictionary& dict
            );

            //- Attempt to read an IOField<Type> and store on objectRegistry
            template<class Type>
            static bool constructIOField
            (
                const word& name,
                token& tok,
                Istream& is,
                objectRegistry& obr
            );

public:

    //- Runtime type information
    TypeName("mappedPatchBase");


    // Constructors

        //- Construct from patch
        explicit mappedPatchBase(const polyPatch&);

        //- Construct with offsetMode=non-uniform
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vectorField& offsets
        );

        //- Construct from offsetMode=uniform
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vector& uniformOffset
        );

        //- Construct from offsetMode=normal and distance
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar normalDistance
        );

        //- Construct from dictionary
        mappedPatchBase(const polyPatch&, const dictionary&);

        //- Construct from dictionary and (collocated) sample mode
        //  (only for nearestPatchFace, nearestPatchFaceAMI, nearestPatchPoint)
        //  Assumes zero offset.
        mappedPatchBase(const polyPatch&, const sampleMode, const dictionary&);

        //- Construct as copy, resetting patch
        mappedPatchBase(const polyPatch&, const mappedPatchBase&);

        //- Construct as copy, resetting patch, map original data
        mappedPatchBase
        (
            const polyPatch&,
            const mappedPatchBase&,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~mappedPatchBase();


    // Member Functions

        // Edit

        void clearOut();

        //- Change to normal offset with given distance
        void setOffset(const scalar normalDist);

        //- Change to uniform offset with value
        void setOffset(const vector& uniformOffset);

        //- Change to non-uniform offsets
        void setOffset(const vectorField& offsets);


        // Access

            //- What to sample
            inline sampleMode mode() const noexcept;

            //- World to sample
            inline const word& sampleWorld() const noexcept;

            //- Region to sample
            inline const word& sampleRegion() const;

            //- Patch (only if NEARESTPATCHFACE)
            inline const word& samplePatch() const;

            //- PatchGroup (only if NEARESTPATCHFACE)
            inline const word& coupleGroup() const;

            //- Return size of mapped mesh/patch/boundary
            inline label sampleSize() const;

            //- Offset vector (from patch faces to destination mesh objects)
            inline const vector& offset() const noexcept;

            //- Offset vectors (from patch faces to destination mesh objects)
            inline const vectorField& offsets() const noexcept;

            //- Get the communicator (worldComm or world-to-world)
            inline label getCommunicator() const;

            //- Identical to getCommunicator()
            inline label comm() const;

            //- Is sample world the local world?
            inline bool sameWorld() const;

            //- Is my world ordered before the sampleWorld?
            inline bool masterWorld() const;

            //- Cached sampleRegion != mesh.name()
            inline bool sameRegion() const noexcept;

            //- Local mesh update time
            inline uniformDimensionedScalarField& updateSampleMeshTime() const;

            //- Sample mesh upate time
            inline uniformDimensionedScalarField& updateMeshTime() const;

            bool upToDate() const;

            //- Return reference to the parallel distribution map
            inline const mapDistribute& map() const;

            //- Return reference to the AMI interpolator
            inline const AMIPatchToPatchInterpolation& AMI
            (
                const bool forceUpdate = false
            ) const;

            //- Is it owner
            inline bool owner() const;

            //- Return a pointer to the AMI projection surface
            const autoPtr<Foam::searchableSurface>& surfPtr() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;


        // Helpers

            //- Get the sample points
            tmp<pointField> samplePoints() const;

            //- Get a point on the face given a face decomposition method:
            //  face-centre-tet : face centre. Returns index of face.
            //  face-planes     : face centre. Returns index of face.
            //  face-diagonal   : intersection of ray from cellcentre to
            //                    facecentre with any of the triangles.
            //                    Returns index (0..size-2) of triangle.
            static pointIndexHit facePoint
            (
                const polyMesh&,
                const label facei,
                const polyMesh::cellDecomposition
            );


        // For database storage

            inline const fileName& sampleDatabasePath() const
            {
                return *sampleDatabasePtr_;
            }

            inline bool sampleDatabase() const
            {
                return bool(sampleDatabasePtr_);
            }

            //- Helper: return path to store data to be sent to processor i
            static fileName sendPath(const fileName& root, const label proci);

            virtual fileName sendPath(const label proci) const;

            //- Helper: return path to store data to be received from
            //- processor i
            static fileName receivePath
            (
                const fileName& root,
                const label proci
            );

            virtual fileName receivePath(const label proci) const;

            //- Lookup (sub)objectRegistry from '/' separated path (relative to
            //- objectRegistry). Creates non-existing intermediate ones.
            static const objectRegistry& subRegistry
            (
                const objectRegistry& obr,
                const fileName& path
            );

            //- Store an IOField on the objectRegistry relative to obr
            template<class Type>
            static void storeField
            (
                objectRegistry& obr,
                const word& fieldName,
                const Field<Type>& values
            );

            //- Convert objectRegistry contents into dictionary
            static void writeDict
            (
                const objectRegistry& obr,
                dictionary& dict
            );

            //- (recursively) construct and register IOFields from dictionary
            static void readDict(const dictionary& d, objectRegistry& obr);


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void distribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void distribute(List<Type>& lst, const CombineOp& cop) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            void reverseDistribute(List<Type>& lst) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            void reverseDistribute(List<Type>& lst, const CombineOp& cop) const;


        // I/O

            //- Write as a dictionary
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
