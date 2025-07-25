/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::ccm::reader

Description
    Reads CCM files as written by PROSTAR/STARCCM

    - arbitrary polyhedral meshs
    - writes interfaces (baffles) to constant/polymesh/interfaces
    - experimental handling of monitoring (internal) boundaries
    - does not handle cyclics. Use createPatch to recreate these
    - does patch names only if they are in the problem description
    - reads cell and face solution data
    - reads Lagrangian data

    The Default_Boundary_Region (region 0) is a special region that serves
    two purposes:
    -# it contains all wall and baffle boundaries that have not otherwise
       been assigned.
    -# it holds the outer bounds of flow domains (fluid/porous/solid)

    The CCM node \c Meshes/FaceBasedTopology/Cells/Interfaces holds the mapping
    of the corresponding mesh faces, which can be used to merge these internal
    boundaries.

    If solid cells exist, there are three possible courses of action:
    -# Remove all solid cells for subsequent flow calculations.
       This is the default.
    -# Treat solid cells like fluid cells, but convert the corresponding
       Default_Boundary_Region to Default_Boundary_Solid for easier
       identification.
       This treatment is useful for visualization purposes.
    -# Move solid cells to a separate mesh region.
       This would be useful for conjugate heat transfer, but
       is not implemented.

    \par Files

    The <tt>constant/remapping</tt> file is an \c IOdictionary that is
    \c READ_IF_PRESENT and can be used to remap certain information. eg,

    \verbatim
        // rename/combine cellTable entries
        //   newName ( listOldNames );
        cellTable
        {
            fluid ( inletRegion outletRegion );
            cat1  ( CAT1 "cat1_(Back|Front|Gamma)" );
        }

        // rename boundary regions
        //   newName oldName;
        boundaryRegion
        {
            inlet_4  inlet_1;
            inlet_5  inlet_2;
            inlet_6  inlet_3;
        }
    \endverbatim

    The <tt>constant/boundaryRegion</tt> file is an \c IOMap<dictionary>
    that is written. It contains the boundary type and names. eg,

    \verbatim
        (
            0
            {
                BoundaryType    wall;
                Label           Default_Boundary_Region;
            }
            1
            {
                BoundaryType    inlet;
                Label           inlet_1;
            }
            ...

            4
            {
                BoundaryType    pressure;
                Label           outlet;
            }
        )
    \endverbatim

    The <tt>constant/cellTable</tt> file is an \c IOMap<dictionary> that is
    written. It contains the cellTable information.
    eg,

    \verbatim
        (
            1
            {
                Label           inletRegion;
                MaterialType    fluid;
                MaterialId      1;
            }
            2
            {
                Label           cat1;
                MaterialType    fluid;
                MaterialId      1;
                PorosityId      1;
            }
            3
            {
                Label           outletRegion;
                MaterialType    fluid;
                MaterialId      1;
            }
        )
    \endverbatim

Note
    this class is still under development
    - any/all of the class names and members may change

SourceFiles
    ccmReader.C
    ccmReaderAux.C
    ccmReaderMesh.C
    ccmReaderSolution.C

\*---------------------------------------------------------------------------*/

#ifndef ccmReader_H
#define ccmReader_H

#include "ccmBase.H"
#include "STARCDCore.H"

#include "labelList.H"
#include "ListOps.H"
#include "polyMesh.H"
#include "boundaryRegion.H"
#include "cellTable.H"
#include "ccmInterfaceDefinitions.H"
#include "ccmSolutionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class ccmID;
class ccmNODE;
class ccmGlobalState;
class reader;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                         Class ccm::reader Declaration
\*---------------------------------------------------------------------------*/

class reader
:
    public base,
    protected fileFormats::STARCDCore
{
public:

    // Forward declarations
    class options;


private:

    // Static Data

        //- cellTable integer options
        static const char* cellTableOpti[];

        //- cellTable string options
        static const char* cellTableOptstr[];


    // Private Data

        //- Reader options
        const std::unique_ptr<options> options_;


        //- Enumeration defining the status of a ccmio node
        enum nodeStatus
        {
            UNKNOWN, BAD, OKAY, READ
        };


        //- status of the geometry (topology) information
        nodeStatus geometryStatus_;

        //- status of the solution (field) information
        nodeStatus solutionStatus_;

        // Persistent data

        //- ccm node: /InterfaceDefinitions
        interfaceDefinitions interfaceDefinitions_;

        //- ccm node: ProblemDescriptions/boundaryRegion
        boundaryRegion boundaryRegion_;

        //- ccm node: ProblemDescriptions/cellType
        cellTable cellTable_;

        //- Number of points
        label nPoints_;

        //- Number of internal faces
        label nInternalFaces_;

        //- Number of faces
        label nFaces_;

        //- Number of cells
        label nCells_;

        //- Minimum mesh data

        //- Points
        pointField points_;

        //- Faces
        faceList faces_;

        //- Face-owner cells
        labelList faceOwner_;

        //- Face-neighbour cells
        labelList faceNeighbour_;

        // List of interface pairs (baffles)
        List<labelPair> bafInterfaces_;

        //- List of interface pairs between mesh regions (domains)
        List<labelPair> domInterfaces_;

        //- Face sets for monitoring
        HashTable<labelList> monitoringSets_;


    // Mesh Maps

        //- Map to original face id
        labelList origFaceId_;

        //- Map to original cell id
        labelList origCellId_;

        //- Cell table id for each cell
        labelList cellTableId_;

        //- List of the original ccm boundary region number
        labelList origBndId_;

        //- Patching and region info
        labelList patchSizes_;

        //- Solution information
        solutionTable   solutionTable_;

        //- Field information
        fieldTable      fieldTable_;

        //- Field information
        fieldTable      lagrangianTable_;


    // Private Member Functions

        //- No copy construct
        reader(const reader&) = delete;

        //- No copy assignment
        void operator=(const reader&) = delete;


        //- Calculate patch starts from given index with current patch sizes
        inline labelList patchStartList(label initial) const;

        //- Report mesh sizes to stdout
        void printSizes() const;

        //- Simulate CCMIOGetEntityIndex for Nodes (not ids)
        int ccmGetEntityIndex(ccmNODE node);

        //- Read string option from specified ccm node
        //  return empty string on failure
        std::string ccmReadNodestr(const char* opt, ccmNODE node);

        //- Read string option from specified ccm node
        //  return empty string on failure
        std::string ccmReadOptstr(const char* opt, ccmID node);

        //- Read map data and check error
        void readMap(const ccmID& mapId, labelList& data);

        //- Determine if the geometry looks good
        bool detectGeometry();

        //- Get interfaces, cellTable and boundaryRegion information
        void readProblemDescription(const ccmID& probNode);

        //- Get /InterfaceDefinitions, used by STARCCM to define in-place
        //  interfaces, etc.
        //  Only handle in-place (IN_PLACE) ones at the moment
        void readInterfaceDefinitions();

        //- Get boundaryRegion information
        void readProblemDescription_boundaryRegion(const ccmID& probNode);

        //- Get cellTable information
        void readProblemDescription_cellTable(const ccmID& probNode);

        //- Read the geometry without any sorting or renumbering
        void readMeshTopology(const scalar scaleFactor=1.0);

        //- Read the vertices
        labelList readVertices
        (
            const ccmID& verticesNode,
            const scalar scaleFactor = 1.0
        );

        //- Read the cells
        void readCells(const ccmID& topoNode);

        //- Read the interfaces
        void readInterfaces(const ccmID& cellsNode);

        //- Read the monitoring
        void readMonitoring(const ccmID& topoId);

        //- Move solid faces from Default_Boundary_Region to
        //  Default_Boundary_Solid
        void juggleSolids();

        //- Remove unwanted fluid/porous/solid regions
        void removeUnwanted();

        //- Remove interfaces with face >= nFace
        void validateInterface(List<labelPair>&);

        //- Renumber interface faces
        void renumberInterfaces(const labelUList& oldToNew);

        //- Remove interfaces between domains (fluid/porosity; fluid/solid, etc)
        //  reorganize baffle interfaces into [0-N/2; N/2-N] lists at the
        //  beginning of the corresponding patch
        void cleanupInterfaces();

        //- Merge the points and faces of in-place conformal interfaces
        void mergeInplaceInterfaces();

        //- Re-order mesh and ensure upper-triangular order
        void reorderMesh();

        //- Write interface (baffle) mapping
        void writeInterfaces(const objectRegistry&) const;

        //- Write List<label> in constant/polyMesh
        void writeMeshLabelList
        (
            const objectRegistry& registry,
            const word& propertyName,
            const labelList& list,
            IOstreamOption streamOpt
        ) const;

        // polyMesh Friend Functions
        void addPatches(polyMesh& mesh) const;

        //- Add faceZones based on monitoring boundary conditions
        void addFaceZones(polyMesh& mesh) const;

        //- Get information about all available solutions
        bool detectSolution();

        //- Get information about available fields
        //  assume that all fields are available for all solution intervals
        void determineFieldInfo
        (
            const ccmID& fieldSetNode,
            fieldTable& table
        );


    // Static Members

        //- Get map of porous regions
        static Map<word> selectPorous(const Map<dictionary>& table);


public:

    // Static Members

        //- Warn about repeated name
        static void warnDuplicates(const word& context, const wordList& lst);


    // Constructors

        //- Open a file for reading
        reader(const fileName& file, const reader::options& opts);


    //- Destructor (closes file)
    ~reader();


    // Member Functions

      // Access

        //- Reference to the reader options
        const reader::options& option() const;


        //- Construct the polyMesh from the read geometry
        //  provide optional remapping dictionary
        autoPtr<polyMesh> mesh
        (
            const objectRegistry& registry,
            const fileName& remappingDictName = fileName::null
        );


        // label nPoints() const { return nPoints_; }
        // label nInternalFaces() const { return nInternalFaces_; }
        // label nFaces() const { return nFaces_; }
        // label nCells() const { return nCells_; }

        // const faceList&   faces() const { return faces_; }
        // const labelList&  faceOwner() const { return faceOwner_; }
        // const labelList&  faceNeighbour() const { return faceNeighbour_; }
        // const pointField& points() const { return points_; }


      // Check

        //- Return true if file has geometry associated with it
        bool hasGeometry();

        //- Return true if file has solutions associated with it
        bool hasSolution();


      // Edit

        //- Remap cellTable and boundaryRegion according to dictionary
        bool remapMeshInfo
        (
            const objectRegistry& registry,
            const fileName& remappingDictName = fileName::null
        );


      // Write

        //- Write the polyMesh
        void writeMesh
        (
            const polyMesh& mesh,
            IOstreamOption streamOpt = IOstreamOption(IOstreamOption::BINARY)
        ) const;

        //- Write cellTable, boundaryRegion and interface information
        void writeAux(const objectRegistry& registry) const;

        //- Detect and read geometry if possible
        bool readGeometry(const scalar scaleFactor = 1.0);

        //- Print general information about the mesh
        void printInfo() const;

        //- Clear out some information after obtaining a polyMesh
        void clearGeom();

        //- Map to original cell Id
        const labelList& origCellId() const
        {
            return origCellId_;
        }

        //- Map to original face Id
        const labelList& origFaceId() const
        {
            return origFaceId_;
        }

        //- Return interface definitions map
        const interfaceDefinitions& interfaceDefinitionsInfo() const
        {
            return interfaceDefinitions_;
        }

        //- Return boundaryRegion table
        const boundaryRegion& boundaryTableInfo() const
        {
            return boundaryRegion_;
        }

        //- Return cell table
        const cellTable& cellTableInfo() const
        {
            return cellTable_;
        }

        //- Return a list of names corresponding to fluids
        Map<word> fluids() const
        {
            return cellTable_.fluids();
        }

        //- Return a list of names corresponding to solids
        Map<word> solids() const
        {
            return cellTable_.solids();
        }

        //- Return table of available solutions
        const solutionTable& solutions()
        {
            detectSolution();
            return solutionTable_;
        }

        //- Return table of available fields
        const fieldTable& fields()
        {
            detectSolution();
            return fieldTable_;
        }

        //- Return table of available lagrangian fields
        const fieldTable& lagrangian()
        {
            detectSolution();
            return lagrangianTable_;
        }

        //- Read solution and field combination
        tmp<scalarField> readField
        (
            const word& solutionName,
            const word& fieldName,
            const bool wallData = false
        );

};


/*---------------------------------------------------------------------------*\
                         Class ccm::reader::options Declaration
\*---------------------------------------------------------------------------*/

class reader::options
{
    // Private data

        //- Keep fluid regions (default true)
        bool keepFluid_;

        //- Keep porous regions (default true)
        bool keepPorous_;

        //- Keep solid regions (default true)
        bool keepSolid_;

        //- Merge in-place interfaces (default true)
        bool mergeInterfaces_;

        //- Rename interface boundaries as InterfaceN_0, InterfaceN_1
        //  (default true)
        bool renameInterfaces_;

        //- Remove baffles by merging their respective faces (default false)
        bool removeBaffles_;

        //- Use numbered names (eg, patch_0, zone_0) instead of human-readable
        //  names (default false)
        bool useNumberedNames_;

        //- merge tolerance for points (default 0.05e-3)
        scalar mergeTol_;

        //- Value to assign for undefined solutions (default: NaN)
        scalar undefScalar_;


public:

    // Constructors

        //- Construct with the defaults
        options();


    // Member Functions

      // Access

        //- Keep fluid regions (default true)
        bool keepFluid() const;

        //- Keep porous regions (default true)
        bool keepPorous() const;

        //- Keep solid regions (default true)
        bool keepSolid() const;

        //- Some region (fluid, porous, solid) is kept.
        bool keptSomeRegion() const;

        //- Merge in-place interfaces (default true)
        bool mergeInterfaces() const;

        //- Rename interface boundaries as InterfaceN_0, InterfaceN_1
        //  (default true)
        bool renameInterfaces() const;

        //- Remove baffles by merging their respective faces (default false)
        bool removeBaffles() const;

        //- Use numbered names (eg, patch_0, zone_0) instead of human-readable
        //  names (default false)
        bool useNumberedNames() const;

        //- Merge tolerance for points (default 0.05e-3)
        scalar mergeTol() const;

        //- Value to assign for undefined solutions (default: NaN)
        scalar undefScalar() const;


      // Edit

        //- Keep fluid regions
        void keepFluid(bool b);

        //- Keep porous regions
        void keepPorous(bool b);

        //- Keep solid regions
        void keepSolid(bool b);

        //- Merge in-place interfaces
        void mergeInterfaces(bool b);

        //- Rename interface boundaries as InterfaceN_0, InterfaceN_1
        void renameInterfaces(bool b);


        //- Remove baffles by merging their respective faces
        void removeBaffles(bool b);

        //- Use numbered names (eg, patch_0, zone_0) instead of human-readable
        void useNumberedNames(bool b);

        //- Merge tolerance for points (default 0.05e-3)
        void mergeTol(const scalar& tol);

        //- Value to assign for undefined solutions (default: NaN)
        void undefScalar(const scalar& val);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
