/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::ccm::writer

Description
    Write OpenFOAM meshes and/or results to CCM format

    - partial support for interfaces
    - no monitoring (internal) boundaries
    - does not handle Lagrangian data

    \par Files

    The <tt>constant/boundaryRegion</tt> and <tt>constant/cellTable</tt> files,
    which are described in ccmReader, are used to construct the CCM
    \c ProblemDescription node.

    The <tt>constant/remapping</tt> file is an \c IOdictionary that is
    \c READ_IF_PRESENT and can be used to remap certain information.
    eg,

    \verbatim
        // map OpenFOAM scalar fields to CCM output fields
        fields
        {
            tracer0
            {
                name    CONC_001;
                Label   "tracer0";
            }
            tracer1
            {
                name    CONC_002;
                Label   "tracer1";
            }
            tracer2
            {
                name    CONC_003;
                Label   "tracer2";
            }
            divPhi
            {
                name    CONC_004;
                Label   "divPhi";
            }
            // an example with units:
            p
            {
                name    P;
                Label   "Pressure";
                units   "Pa";
            }
        }
    \endverbatim

Note
    This class is in development
    - any/all of the class names and members may change

SourceFiles
    ccmWriter.C
    ccmWriterMesh.C
    ccmWriterSolution.C

\*---------------------------------------------------------------------------*/

#ifndef ccmWriter_H
#define ccmWriter_H

#include "ccmBase.H"
#include "STARCDCore.H"

#include "fvMesh.H"
#include "boundaryRegion.H"
#include "cellTable.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class ccmID;
class ccmMaps;
class ccmDimension;
class ccmGlobalState;

/*---------------------------------------------------------------------------*\
                         Class ccm::writer Declaration
\*---------------------------------------------------------------------------*/

class writer
:
    public base,
    protected fileFormats::STARCDCore
{
    // Private Data

    // Static Data

        //- The OpenFOAM -> CCM field mappings
        static dictionary defaultNameMapping;


    // Member Data

        //- MapIds for various components (cell, internalFaces, boundaries)
        std::unique_ptr<ccmMaps> maps_;

        //- mesh reference
        const polyMesh& mesh_;

        //- Lookup between cellModel shape and PROSTAR shape
        const Map<label> prostarShapeLookup_;

        //- ccm node: ProblemDescriptions/boundaryRegion
        boundaryRegion boundaryRegion_;

        // ccm node: ProblemDescriptions/cellType
        cellTable cellTable_;


    // Private Member Functions

        //- No copy construct
        writer(const writer&) = delete;

        //- No copy assignment
        void operator=(const writer&) = delete;


        //- create a filled linear map with 'size' from 'start + 1'
        void addLinearMap
        (
            const string& mapName,
            ccmID& mapId,
            label size,
            label start = 0
        ) const;

        void writeBoundaryRegion(const ccmID& probNode) const;

        void writeCellTable(const ccmID& probNode) const;

        //- write out problem description
        void writeProblem(const ccmID& stateNode) const;

        //- Return the PROSTAR face id for the given cell/face
        label prostarCellFaceId(const label& cellId, const label& faceI) const;

        // write the faces, the number of vertices appears before each entry
        void writeFaces
        (
            const ccmID& nodeId,
            const ccmID& mapId,
            bool  isBoundary,    // boundary or internal faces
            label size,
            label start = 0
        ) const;

        void writeVertices(const ccmID& verticesNode) const;

        // - write internal faces with owner/neighbour
        void writeInternalFaces(const ccmID& topoNode) const;

        // - write boundary faces with owner
        void writeBoundaryFaces(const ccmID& topoNode) const;

        void writeCells(const ccmID& topoNode);

        void writeInterfaces(const ccmID& cellsNode) const;

        bool newFieldNode
        (
            const ccmID& phaseNode,
            const word& fieldName,
            const dictionary& nameMapping,
            const ccmDimension& ccmDim,
            ccmID& fieldNode
        ) const;


        //- Return patch named 'Default_Boundary_Region' or -1 on error
        label findDefaultBoundary() const;

public:

    // Static data members

        //- The name for the topology file reference
        static string defaultMeshName;


    // Constructors

        //- Open a file for writing, with backup/overwrite existing file
        writer
        (
            const fileName& file,
            const polyMesh& mesh,
            const bool backup=true
        );


    //- Destructor (closes file)
    ~writer();


    // Member Functions

      // Write

        //- Write the mesh
        void writeGeometry();

        //- Write the solutions
        //  provide optional remapping dictionary
        void writeSolution
        (
            const IOobjectList& objects,
            const fileName& remappingDictName = fileName::null
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
