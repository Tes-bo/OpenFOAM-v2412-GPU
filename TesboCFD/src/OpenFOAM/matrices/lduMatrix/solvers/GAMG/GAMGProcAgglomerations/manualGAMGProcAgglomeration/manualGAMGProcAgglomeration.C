/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "manualGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        manualGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::manualGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    procAgglomMaps_(controlDict.lookup("processorAgglomeration"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::~manualGAMGProcAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manualGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        forAll(procAgglomMaps_, i)
        {
            const label fineLevelIndex = procAgglomMaps_[i].first();

            if (fineLevelIndex >= agglom_.size())
            {
                WarningInFunction
                    << "Ignoring specification for level " << fineLevelIndex
                    << " since outside agglomeration." << endl;

                continue;
            }

            if (agglom_.hasMeshLevel(fineLevelIndex))
            {
                // Get the fine mesh
                const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
                label nProcs = UPstream::nProcs(levelMesh.comm());

                if (nProcs > 1)
                {
                    // My processor id
                    const label myProcID = Pstream::myProcNo(levelMesh.comm());

                    const List<labelList>& clusters =
                        procAgglomMaps_[i].second();

                    // Coarse to fine master processor
                    labelList coarseToMaster(clusters.size());

                    // Fine to coarse map
                    labelList procAgglomMap(nProcs, -1);

                    // Cluster for my processor (with master index first)
                    labelList agglomProcIDs;



                    forAll(clusters, coarseI)
                    {
                        const labelList& cluster = clusters[coarseI];
                        coarseToMaster[coarseI] = cluster[0];

                        forAll(cluster, i)
                        {
                            procAgglomMap[cluster[i]] = coarseI;
                        }

                        const label masterIndex =
                            cluster.find(coarseToMaster[coarseI]);

                        if (masterIndex == -1)
                        {
                            FatalErrorInFunction
                                << "At level " << fineLevelIndex
                                << " the master processor "
                                << coarseToMaster[coarseI]
                                << " is not in the cluster "
                                << cluster
                                << exit(FatalError);
                        }

                        if (cluster.found(myProcID))
                        {
                            // This is my cluster. Make sure master index is
                            // first
                            agglomProcIDs = cluster;
                            std::swap
                            (
                                agglomProcIDs[0],
                                agglomProcIDs[masterIndex]
                            );
                        }
                    }


                    // Check that we've done all processors
                    if (procAgglomMap.found(-1))
                    {
                        FatalErrorInFunction
                            << "At level " << fineLevelIndex
                            << " processor "
                            << procAgglomMap.find(-1)
                            << " is not in any cluster"
                            << exit(FatalError);
                    }


                    // Communicator for the processor-agglomerated matrix
                    comms_.push_back
                    (
                        UPstream::allocateCommunicator
                        (
                            levelMesh.comm(),
                            coarseToMaster
                        )
                    );

                    // Use processor agglomeration maps to do the actual
                    // collecting
                    if (UPstream::myProcNo(levelMesh.comm()) != -1)
                    {
                        GAMGProcAgglomeration::agglomerate
                        (
                            fineLevelIndex,
                            procAgglomMap,
                            coarseToMaster,
                            agglomProcIDs,
                            comms_.back()
                        );
                    }
                }
            }
        }

        // Print a bit
        if (debug)
        {
            Pout<< nl << "Agglomerated mesh overview" << endl;
            printStats(Pout, agglom_);
        }
    }

    return true;
}


// ************************************************************************* //
