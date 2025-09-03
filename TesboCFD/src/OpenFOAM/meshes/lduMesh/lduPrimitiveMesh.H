/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::lduPrimitiveMesh

Description
    Simplest concrete lduMesh that stores the addressing needed by lduMatrix.

SourceFiles
    lduPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduPrimitiveMesh_H
#define Foam_lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private Data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces_;

        //- Concrete interfaces
        PtrList<const lduInterface> primitiveInterfaces_;

        //- Patch field evaluation schedule
        lduSchedule patchSchedule_;

        //- Communicator to use for any parallel communication
        const label comm_;


    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const PtrList<lduPrimitiveMesh>&);

        label findConnectedInterface
        (
            const lduMesh& myMesh,
            const PtrList<lduPrimitiveMesh>& otherMeshes,
            const labelPairList& procAndInterfaces,

            const label nbrProci,
            const label myRank
        ) const;


        //- No copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&) = delete;

        //- No copy assignment
        void operator=(const lduPrimitiveMesh&) = delete;


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("lduPrimitiveMesh");


    // Constructors

        //- Construct from number of cells
        explicit lduPrimitiveMesh(const label nCells);

        //- Construct from components but without interfaces. Add interfaces
        //  separately using addInterfaces
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            const label comm,
            bool reuse
        );


        //- Add interfaces to a mesh constructed without
        void addInterfaces
        (
            lduInterfacePtrsList& interfaces,
            const lduSchedule& ps
        );

        //- Construct from components and re-use storage.
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            PtrList<const lduInterface>& primitiveInterfaces,
            const lduSchedule& ps,
            const label comm
        );

        //- Construct by combining multiple meshes. The meshes come from
        //  processors procIDs:
        //  procIDs[0] : local processor (myMesh)
        //  procIDs[i] : processor where otherMeshes[i-1] comes from
        //  procAgglomMap : for every processor which processor it agglomerates
        //                  onto. The new processor numbers are in compact
        //                  numbering (so ranks in communicator comm), i.e.
        //                  similar to cell-restrict-addressing.
        //                  We need this information to be able to map
        //                  inter-processor interfaces
        // cellOffsets : for every processor the offset it gets in the mesh
        // faceMap     : for every processor, for every face, the destination
        //               face. Negative for flipped faces.
        // boundaryMap : for every processor, for every patch, -1 or the new
        //               patch index in the mesh.
        // boundaryFaceMap : for every processor, for every patch, for every
        //                   patch face:
        //                   - the new internal face (if boundaryMap=-1)
        //                   - the new patch face (if boundaryMap>=0)
        //                   Faces becoming internal are negative for flipped
        //                   faces.
        lduPrimitiveMesh
        (
            const label comm,
            const labelList& procAgglomMap,

            const labelList& procIDs,
            const lduMesh& myMesh,
            const PtrList<lduPrimitiveMesh>& otherMeshes,

            labelList& cellOffsets,
            labelList& faceOffsets,
            labelListList& faceMap,
            labelListList& boundaryMap,
            labelListListList& boundaryFaceMap
        );


    //- Destructor
    virtual ~lduPrimitiveMesh() = default;


    // Member Functions

        // Access

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const
            {
                return false;
            }

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            //- Return non-const ldu addressing
            virtual lduAddressing& lduAddr()
            {
                return *this;
            }

            //- Return a list of pointers for each patch
            //- with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return a non-const list of pointers for each patch
            //- with only those pointing to interfaces being set
            virtual lduInterfacePtrsList& interfaces()
            {
                return interfaces_;
            }

            //- Return a non-const list of primitive interfaces
            PtrList<const lduInterface>& primitiveInterfaces()
            {
                return primitiveInterfaces_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            //  (reference to cached interfaces)
            const lduInterfacePtrsList& rawInterfaces() const
            {
                return interfaces_;
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return comm_;
            }

            //- Return Lower addressing
            virtual const labelUList& lowerAddr() const
            {
                return lowerAddr_;
            }

            //- Return Upper addressing
            virtual const labelUList& upperAddr() const
            {
                return upperAddr_;
            }

            //- Return non-const Lower addressing
            virtual labelList& lowerAddr()
            {
                return lowerAddr_;
            }

            //- Return non-const Upper addressing
            virtual labelList& upperAddr()
            {
                return upperAddr_;
            }

            //- Return patch addressing
            virtual const labelUList& patchAddr(const label i) const
            {
                return interfaces_[i].faceCells();
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }


        // Helper

            //- Select either mesh0 (meshI is 0) or otherMeshes[meshI-1]
            static const lduMesh& mesh
            (
                const lduMesh& mesh0,
                const PtrList<lduPrimitiveMesh>& otherMeshes,
                const label meshI
            );

            //- Gather meshes from other processors using agglomComm.
            //  Received meshes get GAMGInterface.
            static void gather
            (
                const label agglomComm,
                const lduMesh& mesh,
                PtrList<lduPrimitiveMesh>& otherMeshes
            );

            //- Get non-scheduled send/receive schedule
            template<class ProcPatch>
            static lduSchedule nonBlockingSchedule(const lduInterfacePtrsList&);

            //- Calculate upper-triangular order
            static labelList upperTriOrder
            (
                const label nCells,
                const labelUList& lower,
                const labelUList& upper
            );

            //- Check if in upper-triangular ordering
            static void checkUpperTriangular
            (
                const label size,
                const labelUList& l,
                const labelUList& u
            );

            //- Calculate global cell-cells
            static labelListList globalCellCells
            (
                const lduMesh& mesh,
                const globalIndex& globalNumbering
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduPrimitiveMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
