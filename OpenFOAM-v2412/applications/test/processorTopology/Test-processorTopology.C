/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Test/output processor topology

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "OFstream.H"

#include "openfoam_mpi.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();
    argList::addVerboseOption("Set UPstream::debug level");
    argList::addBoolOption("comm-graph", "Test simple graph communicator");
    argList::addNote
    (
        "Create graph of OpenFOAM mesh connections"
    );

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Only meaningful in parallel"
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createPolyMesh.H"

    // Adjacency table
    const labelListList& connectivity =
        mesh.globalData().topology().procAdjacency();

    if (Pstream::master())
    {
        OFstream os("processorTopology.dot");

        os << "// processorTopology" << nl << nl;
        os.beginBlock("graph");

        forAll(connectivity, proci)
        {
            label nconn = 0;

            for (const label neighProci : connectivity[proci])
            {
                if (proci < neighProci)
                {
                    if (nconn++)
                    {
                        os << "  ";
                    }
                    else
                    {
                        os << indent;
                    }
                    os << proci << " -- " << neighProci;
                }
            }

            if (nconn)
            {
                os  << nl;
            }
        }

        os.endBlock();

        Info<< "Wrote processorTopology graph: "
            << runTime.relativePath(os.name()) << nl;

        Info<< nl
            << "Use neato, circo or fdp graphviz tools" << nl;
    }

    if (Pstream::parRun() && args.found("comm-graph"))
    {
        Info<< nl;

        // Local neighbours
        const labelList& neighbours =
            mesh.globalData().topology().procNeighbours();

        Pout<< "Neigbours: " << flatOutput(neighbours) << endl;

        // As integers values
        List<int> connected(neighbours.size());
        List<int> weights(neighbours.size());
        forAll(neighbours, i)
        {
            connected[i] = neighbours[i];
            weights[i] = 1;
        }

        MPI_Comm topoComm;

        int mpiErrorCode =
        MPI_Dist_graph_create_adjacent
        (
            MPI_COMM_WORLD,
            // Connections into this rank
            connected.size(), connected.cdata(), MPI_UNWEIGHTED,
            // Connections out of this rank
            connected.size(), connected.cdata(), MPI_UNWEIGHTED,
            MPI_INFO_NULL,
            0,   // no reordering (apparently broken anyhow)
            &topoComm
        );

        if (mpiErrorCode)
        {
            FatalError
                << "Failed to create topo communicator. Error:"
                << mpiErrorCode << exit(FatalError);
        }

        int topo_rank = 0;
        int topo_nprocs = 0;
        int topo_inCount = 0;
        int topo_outCount = 0;
        int topo_isWeighted = 0;
        MPI_Comm_rank(topoComm, &topo_rank);
        MPI_Comm_size(topoComm, &topo_nprocs);

        {
            int topo_type = 0;
            MPI_Topo_test(topoComm, &topo_type);

            if (MPI_CART == topo_type)
            {
                Info<< "MPI topology : Cartesian" << endl;
            }
            else if (MPI_GRAPH == topo_type)
            {
                Info<< "MPI topology : Graph" << endl;
            }
            else if (MPI_DIST_GRAPH == topo_type)
            {
                Info<< "MPI topology : Distributed graph" << endl;
            }
            else
            {
                Info<< "MPI topology : None" << endl;
            }
        }

        MPI_Dist_graph_neighbors_count
        (
            topoComm,
            &topo_inCount,
            &topo_outCount,
            &topo_isWeighted
        );

        Pout<< "Topo comm with "
            << topo_rank << " / " << topo_nprocs
            << " from " << connected.size() << flatOutput(connected)
            << " numNbr:" << topo_inCount
            << nl;


        List<int> myPatchIds(neighbours.size());
        forAll(myPatchIds, i)
        {
            // Patches to neighbours
            myPatchIds[i] =
                mesh.globalData().topology().procPatchLookup(neighbours[i]);
        }

        List<int> nbrPatchIds(neighbours.size(), Zero);

        mpiErrorCode = MPI_Neighbor_alltoall
        (
            myPatchIds.data(),
            1,   // one element per neighbour
            MPI_INT,
            nbrPatchIds.data(),
            1,   // one element per neighbour
            MPI_INT,
            topoComm
        );

        if (mpiErrorCode)
        {
            FatalError
                << "MPI Error: " << mpiErrorCode << exit(FatalError);
        }

        Pout<< "proc neighbours:" << flatOutput(neighbours)
            << " my patches:" << flatOutput(myPatchIds)
            << " their patches:" << flatOutput(nbrPatchIds)
            << endl;

        MPI_Comm_free(&topoComm);
    }

    Info<< nl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
