/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::zoneDistribute

Description
    Class for parallel communication in a narrow band. It either provides a Map
    with the neighbouring values of the selected region or returns a Map of the
    required values in global addressing. Also holds a reference to the stencil
    Before the data transfer the communication has to be set up:
    exchangeFields_.setUpCommforZone(interfaceCell_);
    Is used in the plicRDF

    Original code supplied by Henning Scheufler, DLR (2019)

    Additional optimization of processor communication
    provided by Tetsuo AOYAGI, RIST (2022), to use a more compact
    exchange of sizes with an updated version of PstreamBuffers.
    This optimization uses additional sendTo/recvFrom member data
    to track the topological connectivity, acting like an on-the-fly
    sub-communicator, and respects corner connectivity.

    -# Initially topological connections are empty (or all false).
    -# Scan the stencil global cellIds (active zones only) and split
       into sub-lists according the originating processor (the sender).
    -# If an originating processor appears/disappears, need to update
       the connectivity information (requires an all-to-all).
    -# When possible, the topological send/recv is used in PstreamBuffers
       finishedSends (minimizes communication).
    .

SourceFiles
    zoneDistributeI.H
    zoneDistribute.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneDistribute_H
#define Foam_zoneDistribute_H

#include "fvMesh.H"
#include "globalIndex.H"
#include "volFields.H"

#include "zoneCPCStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class zoneDistribute Declaration
\*---------------------------------------------------------------------------*/

class zoneDistribute
:
    public MeshObject<fvMesh, TopologicalMeshObject, zoneDistribute>
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            zoneDistribute
        > MeshObject_type;


    // Private Data

        //- Reference to the zone stencil
        zoneCPCStencil& stencil_;

        //- Global number into index of cells/faces
        const globalIndex& globalNumbering_;

        //- Global cell/face index to send for processor-to-processor comms
        List<labelList> send_;

        //- Parallel [cache]: send connectivity (true/false)
        bitSet sendConnections_;

        //- Parallel [cache]: send data to these ranks
        DynamicList<label> sendProcs_;

        //- Parallel [cache]: recv data from these ranks
        DynamicList<label> recvProcs_;

        //- Persistent set of exchange buffers
        PstreamBuffers pBufs_;


    // Private Member Functions

        //- Return local volField value at (cell or face) index
        template<typename Type>
        Type getLocalValue
        (
            const VolumeField<Type>& phi,
            const label localIdx
        ) const;

        //- Gives patchNumber and patchFaceNumber for a given
        //- Geometric volume field
        template<typename Type>
        Type faceValue
        (
            const VolumeField<Type>& phi,
            const label localIdx
        ) const;


public:

    //- Runtime information
    TypeName("zoneDistribute");


    // Constructors

        //- Construct from fvMesh
        explicit zoneDistribute(const fvMesh&);

        //- Selector
        static zoneDistribute& New(const fvMesh&);


    //- Destructor
    virtual ~zoneDistribute() = default;


    // Member Functions

        //- Update stencil with boolList the size has to match mesh nCells
        void setUpCommforZone(const boolList& zone, bool updateStencil=true);

        //- Updates stencil with boolList the size has to match mesh nCells
        void updateStencil(const boolList& zone);

        //- Stencil reference
        const labelListList& getStencil() noexcept
        {
            return stencil_;
        }

        //- Addressing reference
        const globalIndex& globalNumbering() const noexcept
        {
            return globalNumbering_;
        }

        //- Gives patchNumber and patchFaceNumber for a given
        //- Geometric volume field
        template<typename Type>
        Type getValue
        (
            const VolumeField<Type>& phi,
            const Map<Type>& valuesFromOtherProc,
            const label gblIdx
        ) const;

        //- Returns stencil and provides a Map with globalNumbering
        template<typename Type>
        Map<Field<Type>> getFields
        (
            const boolList& zone,
            const VolumeField<Type>& phi
        );

        //- Returns stencil and provides a Map with globalNumbering
        template<typename Type>
        Map<Type> getDatafromOtherProc
        (
            const boolList& zone,
            const VolumeField<Type>& phi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "zoneDistributeI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
