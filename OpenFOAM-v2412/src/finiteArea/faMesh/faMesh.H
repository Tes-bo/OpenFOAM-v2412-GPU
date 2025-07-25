/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::faMesh

Description
    Finite area mesh (used for 2-D non-Euclidian finite area method)
    defined using a \em patch of faces on a polyMesh
    (ie, uindirectPrimitivePatch).

    The ordering of faces and points on the faMesh corresponds to
    the localFaces and localPoints as per Foam::PrimitivePatch but
    the edge addressing is handled slightly differently.
    The internal edges of the faMesh will generally correspond identically
    to the internalEdges of the PrimitivePatch (may change in the future)
    but the boundaryEdges will be reordered compared to the PrimitivePatch
    to allow edge boundary slices to be obtained.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C
    faMeshPatches.C
    faMeshTopology.C
    faMeshUpdate.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMesh_H
#define Foam_faMesh_H

#include "polyMesh.H"
#include "lduMesh.H"
#include "faBoundaryMesh.H"
#include "edgeList.H"
#include "faceList.H"
#include "primitiveFieldsFwd.H"
#include "DimensionedField.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "indirectPrimitivePatch.H"
#include "edgeInterpolation.H"
#include "labelIOList.H"
#include "FieldFields.H"
#include "faGlobalMeshData.H"
#include "faSchemes.H"
#include "faSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMeshBoundaryHalo;
class faMeshLduAddressing;
class faMeshMapper;
class faPatchData;

/*---------------------------------------------------------------------------*\
                       Class faMeshRegistry Declaration
\*---------------------------------------------------------------------------*/

//- The objectRegistry for faMesh.
//  This is a separate class to ensure it will be fully constructed before
//  other data members use it, which ensures that its virtual methods are
//  callable during construction (gcc).
class faMeshRegistry
:
    public objectRegistry
{
public:

    // Constructors

        //- Construct an objectRegistry for given area region name.
        //  Uses faMeshesRegistry internally
        faMeshRegistry(const word& areaRegion, const polyMesh& mesh);


    // Member Functions

        //- True - thisDb() is a valid registry
        virtual bool hasDb() const { return true; }

        //- Reference to the mesh database
        virtual const objectRegistry& thisDb() const { return *this; }

        //- Local directory path of the objectRegistry relative to Time
        //- with override for the single-region case
        virtual const fileName& dbDir() const;
};


/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh
:
    public faMeshRegistry,
    public lduMesh,
    public faSchemes,
    public faSolution,
    public edgeInterpolation  // may need input from faSchemes
{
    // Private (internal) classes/structures

        //- A (proc, patchi, patchEdgei, meshFacei) tuple used internally
        //- for managing patch/patch bookkeeping during construction.
        //  Finite-area patches are stored with negated indices (offset -2),
        //  which makes them readily identifiable and always sort before normal
        //  patches.
        struct patchTuple
        :
            public FixedList<label, 4>
        {

        // Constructors

            // Inherit constructors
            using FixedList<label, 4>::FixedList;

            //- Default construct as 'invalid'
            patchTuple()
            {
                clear();
            }


        // Static Member Functions

            // Globally consistent ordering
            //
            // 1. sort left/right as lower/higher processor connection
            // 2. sort by proc/patch/patch index
            static void sort(UList<Pair<patchTuple>>& list)
            {
                for (auto& tuples : list)
                {
                    tuples.sort();
                }
                Foam::stableSort(list);
            }


        // Member Functions

            //- Reset to 'invalid'
            void clear()
            {
                procNo(-1);
                patchi(labelMax);
                patchEdgei(-1);
                meshFacei(-1);
            }

            //- Valid if proc and patch (or patch edge) are non-negative
            bool valid() const
            {
                return (procNo() >= 0 && patchi() != -1);
            }

            // Processor is the first sort index
            label procNo() const { return (*this)[0]; }
            void procNo(label val) { (*this)[0] = val; }

            // PatchId is the second sort index (finiteArea patches are < -1)
            label patchi() const { return (*this)[1]; }
            void patchi(label val) { (*this)[1] = val; }

            // The patch edge index (on the finiteArea patch)
            // is the third sort index
            label patchEdgei() const { return (*this)[2]; }
            void patchEdgei(label val) { (*this)[2] = val; }

            // The processor-local mesh face is the fourth sort index
            label meshFacei() const { return (*this)[3]; }
            void meshFacei(label val) { (*this)[3] = val; }

            //- Return the real patchId
            label realPatchi() const
            {
                const label id = patchi();
                return (id < -1 ? -(id + 2) : id);
            }

            //- Set patchId as finiteArea
            void faPatchi(label val)
            {
                patchi(-(val + 2));
            }

            //- Considered to be finiteArea if (patchi < -1)
            bool is_finiteArea() const noexcept
            {
                return (patchi() < -1);
            }

            //- Considered to be processor local
            bool is_localProc() const
            {
                return (procNo() == UPstream::myProcNo());
            }
        };


    // Private Data

        //- Face labels
        labelIOList faceLabels_;

        //- Boundary mesh
        faBoundaryMesh boundary_;


    // Primitive mesh data

        //- Edges, addressing into local point list
        edgeList edges_;

        //- Edge owner
        labelList edgeOwner_;

        //- Edge neighbour
        labelList edgeNeighbour_;


    // Primitive size data

        //- Total number of points
        mutable label nPoints_;

        //- Total number of edges
        mutable label nEdges_;

        //- Number of internal edges
        mutable label nInternalEdges_;

        //- Number of faces
        mutable label nFaces_;


    // Communication support, updating

        //- Communicator used for parallel communication
        label comm_;

        //- Current time index for motion
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;


    // Demand-driven data

        //- Primitive patch
        mutable std::unique_ptr<uindirectPrimitivePatch> patchPtr_;

        //- List of poly patch/patch-face for faceLabels
        mutable std::unique_ptr<List<labelPair>> polyPatchFacesPtr_;

        //- List of polyPatch ids used by the area mesh
        mutable std::unique_ptr<labelList> polyPatchIdsPtr_;

        //- List of proc/mesh-face for boundary edge neighbours
        mutable std::unique_ptr<List<labelPair>> bndConnectPtr_;

        //- Ldu addressing data
        mutable std::unique_ptr<faMeshLduAddressing> lduPtr_;


    // Geometric Data

        //- Face areas
        mutable std::unique_ptr<DimensionedField<scalar, areaMesh>> SPtr_;

        //- Face areas old time level
        mutable std::unique_ptr<DimensionedField<scalar, areaMesh>> S0Ptr_;

        //- Face areas old-old time level
        mutable std::unique_ptr<DimensionedField<scalar, areaMesh>> S00Ptr_;

        //- Patch starts in the edge list
        mutable std::unique_ptr<labelList> patchStartsPtr_;

        //- Edge length vectors
        mutable std::unique_ptr<edgeVectorField> LePtr_;

        //- Mag edge length vectors
        mutable std::unique_ptr<edgeScalarField> magLePtr_;

        //- Face centres
        mutable std::unique_ptr<areaVectorField> faceCentresPtr_;

        //- Edge centres
        mutable std::unique_ptr<edgeVectorField> edgeCentresPtr_;

        //- Face area normals
        mutable std::unique_ptr<areaVectorField> faceAreaNormalsPtr_;

        //- Edge area normals
        mutable std::unique_ptr<edgeVectorField> edgeAreaNormalsPtr_;

        //- Point area normals
        mutable std::unique_ptr<vectorField> pointAreaNormalsPtr_;

        //- Face curvatures
        mutable std::unique_ptr<areaScalarField> faceCurvaturesPtr_;

        //- Edge transformation tensors
        mutable std::unique_ptr<FieldField<Field, tensor>>
            edgeTransformTensorsPtr_;

        //- Whether point normals must be corrected for a patch
        mutable std::unique_ptr<boolList> correctPatchPointNormalsPtr_;


    // Other mesh-related data

        //- Parallel info
        mutable autoPtr<faGlobalMeshData> globalMeshDataPtr_;

        //- Mapping/swapping for boundary edge halo neighbours
        mutable std::unique_ptr<faMeshBoundaryHalo> haloMapPtr_;

        //- Face centres for boundary edge halo neighbours
        mutable std::unique_ptr<pointField> haloFaceCentresPtr_;

        //- Face normals for boundary edge halo neighbours
        mutable std::unique_ptr<vectorField> haloFaceNormalsPtr_;


    // Static Private Data

        //- The prefix to local: %finite-area
        static const word prefix_;

        //- Quadrics fit for pointAreaNormals (experimental)
        static const int quadricsFit_;


    // Private Member Functions

        //- Set indirect patch, removing any old one.
        //  No communication
        void initPatch() const;

        //- Set primitive mesh data.
        //  No communication
        void setPrimitiveMeshData();

        //- Get list of (proc/patchi/patchEdgei/meshFacei) tuple pairs in an
        //- globally consistent ordering
        List<Pair<patchTuple>> getBoundaryEdgeConnections() const;

        //- Determine the boundary edge neighbour connections
        void calcBoundaryConnections() const;

        //- Define boundary edge neighbours (proc/face) based on
        //- gathered topology information
        void setBoundaryConnections
        (
            const List<Pair<patchTuple>>& bndEdgeConnections
        ) const;


        // Private Member Functions to calculate demand driven data

            //- Calculate ldu addressing
            void calcLduAddressing() const;

            //- Calculate patch starts in the edge list
            void calcPatchStarts() const;

            //- Calculate which polyPatches, polyPatch/local-face
            //- are related to the areaMesh
            void calcWhichPatchFaces() const;

            //- Calculate the edge normals (direct calculation)
            //- with flat boundary addressing
            //  Possible communication
            tmp<vectorField> calcRawEdgeNormals(int calcType) const;

            //- Calculate edge lengths
            //  Triggers communication via calcEdgeAreaNormals
            void calcLe() const;

            //- Calculate mag edge lengths
            //  No communication
            void calcMagLe() const;

            //- Calculate face centres
            //  No communication
            void calcFaceCentres() const;

            //- Calculate edge centres
            //  No communication
            void calcEdgeCentres() const;

            //- Calculate face areas
            //  No communication
            void calcS() const;

            //- Calculate face area normals
            //  Triggers communication via calcEdgeAreaNormals
            void calcFaceAreaNormals() const;

            //- Calculate edge area normals
            //  Triggers communication via pointAreaNormals
            void calcEdgeAreaNormals() const;

            //- Calculate point area normals
            //  Triggers communication for processor patches
            void calcPointAreaNormals(vectorField& result) const;

            //- Calculate point area normals by quadrics fit
            void calcPointAreaNormalsByQuadricsFit(vectorField& result) const;

            //- Calculate face curvatures
            //  Triggers communication via edgeLengthCorrection
            void calcFaceCurvatures() const;

            //- Calculate edge transformation tensors
            void calcEdgeTransformTensors() const;

            //- Clear geometry but not the face areas
            void clearGeomNotAreas() const;

            //- Has halo face centers/normals
            bool hasHaloFaceGeometry() const noexcept;

            //- Clear boundary halo information
            void clearHalo() const;

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

            //- Clear demand-driven data
            void clearOut() const;


    // Halo handling

        //- Calculate halo centres/normals
        void calcHaloFaceGeometry() const;


    // Helpers

        //- Create a single patch
        faPatchList createOnePatch
        (
            const word& patchName,
            const word& patchType = ""
        ) const;

        //- Create list of patches from boundary definition
        faPatchList createPatchList
        (
            const dictionary& bndDict,
            const word& emptyPatchName = "",
            const dictionary* defaultPatchDefinition = nullptr
        ) const;


        //- Fatal error if edge labels are out of range
        void checkBoundaryEdgeLabelRange(const labelUList& edgeLabels) const;

        //- Extract list from contiguous (unordered) boundary data
        //- to the locally sorted order.
        template<class T>
        List<T> boundarySubset
        (
            const UList<T>& bndField,
            const labelUList& edgeLabels
        ) const
        {
            #ifdef FULLDEBUG
            checkBoundaryEdgeLabelRange(edgeLabels);
            #endif
            // Like UIndirectList but with an offset
            List<T> result(edgeLabels.size());
            forAll(edgeLabels, i)
            {
                result[i] = bndField[edgeLabels[i] - nInternalEdges_];
            }
            return result;
        }


    // Static Functions

        //- Test if faSchemes/faSolution files are available
        static bool hasSystemFiles
        (
            const word& meshName,
            const polyMesh& pMesh
        );

        //- Test if mesh files needed for read construction are available
        static bool hasMeshFiles
        (
            const word& meshName,
            const polyMesh& pMesh
        );

public:

    // Public Typedefs

        //- The mesh type
        typedef faMesh Mesh;

        //- The boundary type associated with the mesh
        typedef faBoundaryMesh BoundaryMesh;


    // Tuning switches

        //- Geometry treatment
        static int geometryOrder_;


    //- Runtime type information
    TypeName("faMesh");

        //- The prefix to the parent registry name: %finite-area
        static const word& prefix() noexcept;

        //- The mesh sub-directory name (usually "faMesh")
        static word meshSubDir;


    // Generated Methods

        //- No copy construct
        faMesh(const faMesh&) = delete;

        //- No copy assignment
        void operator=(const faMesh&) = delete;


    // Constructors

        //- Construct zero-sized from polyMesh
        //  Boundary is added using addFaPatches() member function
        faMesh(const word& meshName, const polyMesh&, Foam::zero);

        //- Construct zero-sized from polyMesh
        //  Boundary is added using addFaPatches() member function
        faMesh(const polyMesh&, Foam::zero);

        //- Read construct from polyMesh, using its IOobject properties
        faMesh(const word& meshName, const polyMesh&, const bool doInit = true);

        //- Read construct from polyMesh, using its IOobject properties
        explicit faMesh(const polyMesh&, const bool doInit = true);

        //- Construct from components (face labels) without boundary,
        //- using specified read properties.
        //  Boundary is added using addFaPatches() member function.
        faMesh
        (
            const word& meshName,
            const polyMesh& pMesh,
            labelList&& faceLabels,
            IOobjectOption ioOpt
        );

        //- Construct from components (face labels) without boundary,
        //- using IOobject properties from polyMesh.
        //  Boundary is added using addFaPatches() member function.
        faMesh(const word& meshName, const polyMesh&, labelList&& faceLabels);

        //- Construct from components (face labels) without boundary,
        //- using specified read properties.
        //  Boundary is added using addFaPatches() member function.
        faMesh(const polyMesh&, labelList&& faceLabels, IOobjectOption ioOpt);

        //- Construct from components (face labels) without boundary,
        //- using IOobject properties from polyMesh.
        //  Boundary is added using addFaPatches() member function.
        faMesh(const polyMesh&, labelList&& faceLabels);

        //- Construct as copy (for dictionaries) and zero-sized
        //- without boundary.
        //  Boundary is added using addFaPatches() member function
        faMesh(const word& meshName, const faMesh& baseMesh, Foam::zero);

        //- Construct as copy (for dictionaries) and zero-sized
        //- without boundary, using IOobject properties from polyMesh.
        //  Boundary is added using addFaPatches() member function
        faMesh(const faMesh& baseMesh, Foam::zero);

        //- Construct as copy (for dictionaries) and faceLabels
        //- without boundary, using specified read properties.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const word& meshName,
            const faMesh& baseMesh,
            labelList&& faceLabels,
            IOobjectOption ioOpt
        );

        //- Construct as copy (for dictionaries) and faceLabels
        //- without boundary, using specified read properties.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const faMesh& baseMesh,
            labelList&& faceLabels,
            IOobjectOption ioOpt
        );

        //- Construct as copy (for dictionaries) and faceLabels
        //- without boundary, using IOobject properties from polyMesh.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const word& meshName,
            const faMesh& baseMesh,
            labelList&& faceLabels
        );

        //- Construct as copy (for dictionaries) and faceLabels
        //- without boundary, using read properties from baseMesh.
        //  Boundary is added using addFaPatches() member function
        faMesh(const faMesh& baseMesh, labelList&& faceLabels);

        //- Construct from single polyPatch
        faMesh
        (
            const word& meshName,
            const polyPatch& pp,
            const bool doInit = true
        );

        //- Construct from single polyPatch
        explicit faMesh(const polyPatch& pp, const bool doInit = true);

        //- Construct from definition
        faMesh
        (
            const word& meshName,
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition,
            const bool doInit = true
        );

        //- Construct from definition
        faMesh
        (
            const polyMesh& pMesh,
            const dictionary& faMeshDefinition,
            const bool doInit = true
        );


    //- Destructor
    virtual ~faMesh();


    // Static Functions

        //- Return the current geometry treatment
        //  A zero level or negative is with restricted neighbour information
        static int geometryOrder() noexcept
        {
            return geometryOrder_;
        }

        //- Set the preferred geometry treatment
        //  \return the previous value
        static int geometryOrder(int order) noexcept
        {
            int old(geometryOrder_);
            geometryOrder_ = order;
            return old;
        }

        //- Read construction from polyMesh if all files are available
        static autoPtr<faMesh> TryNew
        (
            const word& meshName,
            const polyMesh& pMesh
        );

        //- Read construction from polyMesh if all files are available
        static autoPtr<faMesh> TryNew(const polyMesh& pMesh);


    // Member Functions

    // Topological Change

        //- Add boundary patches. Constructor helper
        void addFaPatches
        (
            faPatchList& plist,
            const bool validBoundary = true
        );

        //- Add boundary patches. Constructor helper
        void addFaPatches
        (
            const List<faPatch*>& p,
            const bool validBoundary = true
        );

        //- Initialise non-demand-driven data etc
        bool init(const bool doInit);

        //- Processor/processor synchronisation for geometry fields.
        //  Largely internal use only (slightly hacky).
        void syncGeom();


    // Database

        //- The parent registry containing all finite-area meshes
        //- on the polyMesh.
        static const objectRegistry* registry(const polyMesh& pMesh);

        //- The single-region finite-area region on the polyMesh.
        //- Uses lookupObject semantics - Fatal if non-existent
        static const faMesh& mesh(const polyMesh& pMesh);

        //- The single-region or specified finite-area region on the polyMesh.
        //- Uses lookupObject semantics - Fatal if non-existent
        static const faMesh& mesh(const polyMesh&, const word& areaRegion);

        //- Return access to polyMesh
        const polyMesh& mesh() const;

        //- Return the local mesh directory (dbDir()/meshSubDir)
        fileName meshDir() const;

        //- Return reference to time
        const Time& time() const;

        //- Return the current instance directory for points
        //  Used in the construction of geometric mesh data dependent
        //  on points
        const fileName& pointsInstance() const;

        //- Return the current instance directory for faces
        const fileName& facesInstance() const;

        //- Const reference to the mesh and solver state data
        const meshState& data() const { return mesh().data(); }

        //- Reference to the mesh and solver state data
        meshState& data() { return const_cast<polyMesh&>(mesh()).data(); }


    // Parallel

        //- Return communicator used for parallel communication
        label comm() const noexcept { return comm_; }

        //- Return communicator used for parallel communication
        label& comm() noexcept { return comm_; }

        //- Is demand-driven parallel info available?
        bool hasGlobalData() const noexcept;

        //- Return parallel info (demand-driven)
        const faGlobalMeshData& globalData() const;


    // Solution Control

        //- Non-null if faSchemes exists (can test as bool).
        const faSchemes* hasSchemes() const;

        //- Non-null if faSolution exists (can test as bool).
        const faSolution* hasSolution() const;

        //- Read-access to the faSchemes controls
        const faSchemes& schemes() const;

        //- Read/write-access to the faSchemes controls
        faSchemes& schemes();

        //- Read-access to the faSolution controls
        const faSolution& solution() const;

        //- Read/write-access to the faSolution controls
        faSolution& solution();


    // Access: Mesh size parameters

        //- Number of local mesh points
        inline label nPoints() const noexcept;

        //- Number of local mesh edges
        inline label nEdges() const noexcept;

        //- Number of internal faces
        inline label nInternalEdges() const noexcept;

        //- Number of boundary edges (== nEdges - nInternalEdges)
        inline label nBoundaryEdges() const noexcept;

        //- Number of patch faces
        inline label nFaces() const noexcept;


    // Access: Primitive mesh data

        //- Return local points
        inline const pointField& points() const;

        //- Return local edges with reordered boundary
        inline const edgeList& edges() const noexcept;

        //- Sub-list of local internal edges
        inline const edgeList::subList internalEdges() const;

        //- Return local faces
        inline const faceList& faces() const;

        //- Edge owner addressing
        inline const labelList& edgeOwner() const noexcept;

        //- Edge neighbour addressing
        inline const labelList& edgeNeighbour() const noexcept;

        //- True if the internalEdges use an ordering that does not
        //- correspond 1-to-1 with the patch internalEdges.
        inline bool hasInternalEdgeLabels() const noexcept;


    // Database

        //- True - thisDb() is a valid registry
        virtual bool hasDb() const { return true; }

        //- Reference to the mesh database
        virtual const objectRegistry& thisDb() const
        {
            return faMeshRegistry::thisDb();
        }

        //- Local directory path of the objectRegistry relative to Time
        //- with override for the single-region case
        virtual const fileName& dbDir() const
        {
            return faMeshRegistry::dbDir();
        }

        //- Name function is needed to disambiguate those inherited
        //- from base classes
        const word& name() const
        {
            return faMeshRegistry::thisDb().name();
        }


    // Regions

        //- Local registry directory path (relative to Time) for specified
        //- area mesh (of a single-region volume mesh)
        static fileName dbDir(const word& areaRegion);

        //- Local registry directory path (relative to Time) for specified
        //- volume mesh and area mesh combination
        static fileName dbDir(const word& volRegion, const word& areaRegion);

        //- Local registry directory path (relative to Time) for specified
        //- volume mesh and area mesh combination
        static fileName dbDir
        (
            const polyMesh& pMesh,
            const word& areaRegion = word::null
        );

        //- The local mesh directory name (eg, "faMesh") for specified
        //- area mesh (of a single-region volume mesh)
        static fileName meshDir(const word& areaRegion);

        //- The local mesh directory name (eg, "faMesh") for specified
        //- volume mesh and area mesh combination
        static fileName meshDir(const word& volRegion, const word& areaRegion);

        //- The local mesh directory name (eg, "faMesh") for specified
        //- volume mesh and area mesh combination
        static fileName meshDir
        (
            const polyMesh& pMesh,
            const word& areaRegion = word::null
        );

        //- The mesh region name or word::null if polyMesh::defaultRegion
        const word& regionName() const;


        // Access

            //- Return constant reference to boundary mesh
            inline const faBoundaryMesh& boundary() const noexcept;

            //- Return the underlying polyMesh face labels
            inline const labelList& faceLabels() const noexcept;

            //- The area-face corresponding to the mesh-face, -1 if not found
            inline label whichFace(const label meshFacei) const;

            //- The polyPatches related to the areaMesh, in sorted order
            inline const labelList& whichPolyPatches() const;

            //- The polyPatch/local-face for each faceLabels()
            inline const List<labelPair>& whichPatchFaces() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Internal face owner
            const labelUList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- True if given edge label is internal to the mesh
            bool isInternalEdge(const label edgeIndex) const noexcept
            {
                return (edgeIndex < nInternalEdges_);
            }

            //- List of proc/face for the boundary edge neighbours
            //- using primitive patch edge numbering.
            inline const List<labelPair>& boundaryConnections() const;

            //- Boundary edge neighbour processors
            //- (does not include own proc)
            labelList boundaryProcs() const;

            //- List of proc/size for the boundary edge neighbour processors
            //- (does not include own proc)
            List<labelPair> boundaryProcSizes() const;

            //- Mapping/swapping for boundary halo neighbours
            const faMeshBoundaryHalo& boundaryHaloMap() const;

            //- Face centres of boundary halo neighbours
            const pointField& haloFaceCentres() const;

            //- Face unit-normals of boundary halo neighbours
            const vectorField& haloFaceNormals() const;

            //- Face centres of boundary halo neighbours for specified patch
            tmp<pointField> haloFaceCentres(const label patchi) const;

            //- Face unit-normals of boundary halo neighbours for specified patch
            tmp<vectorField> haloFaceNormals(const label patchi) const;


    // Interfacing

        //- The volume (owner) cells associated with the area-mesh
        labelList faceCells() const;


    // Storage Management

        //- Remove all files from mesh instance
        void removeFiles(const fileName& instanceDir) const;

        //- Remove all files from mesh instance()
        void removeFiles() const;


        //- Has face areas: S()
        bool hasFaceAreas() const noexcept { return bool(SPtr_); }

        //- Has face centres: areaCentres()
        bool hasAreaCentres() const noexcept { return bool(faceCentresPtr_); }

        //- Has edge centres: edgeCentres()
        bool hasEdgeCentres() const noexcept { return bool(edgeCentresPtr_); }

        //- Has edge length vectors: Le()
        bool hasLe() const noexcept { return bool(LePtr_); }

        //- Has edge length magnitudes: magLe()
        bool hasMagLe() const noexcept { return bool(magLePtr_); }

        //- Has face area normals: faceAreaNormals()
        bool hasFaceAreaNormals() const noexcept
        {
            return bool(faceAreaNormalsPtr_);
        }

        //- Has edge area normals: edgeAreaNormals()
        bool hasEdgeAreaNormals() const noexcept
        {
            return bool(edgeAreaNormalsPtr_);
        }

        //- Has point area normals: pointAreaNormals()
        bool hasPointAreaNormals() const noexcept
        {
            return bool(pointAreaNormalsPtr_);
        }

        //- Has face curvatures: faceCurvatures()
        bool hasFaceCurvatures() const noexcept
        {
            return bool(faceCurvaturesPtr_);
        }

        //- Has patch point normals corrections
        bool hasPatchPointNormalsCorrection() const noexcept
        {
            return bool(correctPatchPointNormalsPtr_);
        }


    // Mesh motion and morphing

        //- Is mesh moving
        bool moving() const
        {
            return mesh().moving();
        }

        //- Update after mesh motion
        virtual bool movePoints();

        //- Update after topo change
        virtual void updateMesh(const mapPolyMesh&);


        // Mapping

            //- Map all fields in time using given map.
            virtual void mapFields(const faMeshMapper& mapper) const;

            //- Map face areas in time using given map.
            virtual void mapOldAreas(const faMeshMapper& mapper) const;


        // Demand-driven data

            //- Return constant reference to primitive patch
            inline const uindirectPrimitivePatch& patch() const;

            //- Return reference to primitive patch
            inline uindirectPrimitivePatch& patch();

            //- Return patch starts
            const labelList& patchStarts() const;

            //- Return edge length vectors
            const edgeVectorField& Le() const;

            //- Return edge length magnitudes
            const edgeScalarField& magLe() const;

            //- Return normalised edge length vectors
            tmp<edgeVectorField> unitLe() const;

            //- Return face centres as areaVectorField
            const areaVectorField& areaCentres() const;

            //- Return edge centres as edgeVectorField
            const edgeVectorField& edgeCentres() const;

            //- Return face areas
            const DimensionedField<scalar, areaMesh>& S() const;

            //- Return old-time face areas
            const DimensionedField<scalar, areaMesh>& S0() const;

            //- Return old-old-time face areas
            const DimensionedField<scalar, areaMesh>& S00() const;

            //- Return face area normals
            const areaVectorField& faceAreaNormals() const;

            //- Return edge area normals
            const edgeVectorField& edgeAreaNormals() const;

            //- Return point area normals
            const vectorField& pointAreaNormals() const;

            //- Return face curvatures
            const areaScalarField& faceCurvatures() const;

            //- Return edge transformation tensors
            const FieldField<Field, tensor>& edgeTransformTensors() const;

            //- Return internal point labels
            labelList internalPoints() const;

            //- Return boundary point labels
            labelList boundaryPoints() const;

            //- Return edge length correction
            tmp<edgeScalarField> edgeLengthCorrection() const;

            //- Whether point normals should be corrected for a patch
            bool correctPatchPointNormals(const label patchID) const;

            //- Set whether point normals should be corrected for a patch
            boolList& correctPatchPointNormals() const;


        // Storage management


        //- Write mesh
        virtual bool write(const bool writeOnProc = true) const;


    // Member Operators

        bool operator!=(const faMesh& m) const;

        bool operator==(const faMesh& m) const;


    // Housekeeping

        //- No call operator. Prior to 2312 was used to obtain polyMesh
        void operator()() const = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faMeshI.H"

#ifdef NoRepository
    #include "faPatchFaMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
