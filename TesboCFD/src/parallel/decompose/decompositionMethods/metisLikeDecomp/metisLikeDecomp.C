/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "metisLikeDecomp.H"
#include "Time.H"
#include "globalIndex.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisLikeDecomp::decomposeGeneral
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    if (!UPstream::parRun())
    {
        // Treat zero-sized weights as uniform weighting
        return decomposeSerial
        (
            adjncy,
            xadj,
            cWeights,
            decomp
        );
    }

    if (debug)
    {
        Info<< type() << "Decomp : running in parallel."
            << " Decomposing all of graph on master processor." << endl;
    }

    // Protect against zero-sized offset list
    const label numCells = Foam::max(0, (xadj.size()-1));

    const globalIndex globalAdjncy(adjncy.size());
    const globalIndex globalCells(numCells);

    List<label> allAdjncy(globalAdjncy.gather(adjncy));

    // Gathering xadj to master is similar to globalIndex gather()
    // except for the following:
    //
    //   - gathered list is size+1
    //   - apply local to global renumbering

    const UPstream::commsTypes commsType = UPstream::commsTypes::nonBlocking;
    const label startOfRequests = UPstream::nRequests();


    List<label> allXadj;
    if (UPstream::master())
    {
        allXadj.resize(globalCells.totalSize()+1);
        allXadj.back() = globalAdjncy.totalSize();  // Final end offset

        // My values - no renumbering required
        SubList<label>(allXadj, globalCells.localSize(0)) =
            SubList<label>(xadj, globalCells.localSize(0));

        for (const int proci : globalCells.subProcs())
        {
            SubList<label> procSlot(allXadj, globalCells.range(proci));

            if (procSlot.empty())
            {
                // Nothing to do
            }
            else
            {
                UIPstream::read
                (
                    commsType,
                    proci,
                    procSlot.data_bytes(),
                    procSlot.size_bytes(),
                    UPstream::msgType(),
                    UPstream::worldComm
                );
            }
        }
    }
    else
    {
        // Send my part of the graph (local numbering)

        if (!numCells)
        {
            // Nothing to do
        }
        else
        {
            SubList<label> procSlot(xadj, numCells);

            UOPstream::write
            (
                commsType,
                UPstream::masterNo(),
                procSlot.cdata_bytes(),
                procSlot.size_bytes(),
                UPstream::msgType(),
                UPstream::worldComm
            );
        }
    }

    // Wait for outstanding requests
    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        UPstream::waitRequests(startOfRequests);
    }

    // Local to global renumbering
    if (UPstream::master())
    {
        for (const int proci : globalCells.subProcs())
        {
            SubList<label> procSlot(allXadj, globalCells.range(proci));

            globalAdjncy.inplaceToGlobal(proci, procSlot);
        }
    }

    // Uniform weighting if weights are empty or poorly sized

    List<scalar> allWeights;
    if (returnReduceAnd(cWeights.size() == globalCells.localSize()))
    {
        // ie, hasWeights
        allWeights = globalCells.gather(cWeights);
    }


    // Global decomposition
    labelList allDecomp;

    if (UPstream::master())
    {
        decomposeSerial
        (
            allAdjncy,
            allXadj,
            allWeights,
            allDecomp
        );

        allAdjncy.clear();    // Not needed anymore
        allXadj.clear();      // ...
        allWeights.clear();   // ...
    }

    // The processor-local decomposition (output)
    decomp.resize_nocopy(globalCells.localSize());
    globalCells.scatter(allDecomp, decomp);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisLikeDecomp::metisLikeDecomp(const label numDomains)
:
    decompositionMethod(numDomains),
    coeffsDict_(dictionary::null)
{}


Foam::metisLikeDecomp::metisLikeDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    if (!points.empty() && (points.size() != mesh.nCells()))
    {
        FatalErrorInFunction
            << "Number of cell centres (" << points.size()
            << ") != number of cells (" << mesh.nCells() << ")"
            << exit(FatalError);
    }

    // Global mesh connectivity
    CompactListList<label> cellCells;
    globalMeshData::calcCellCells(mesh, cellCells, true);

    labelList decomp;
    decomposeGeneral
    (
        cellCells.values(),
        cellCells.offsets(),
        pointWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
) const
{
    if (agglom.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Agglomeration size (" << agglom.size()
            << ") != number of cells (" << mesh.nCells() << ")"
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells;
    globalMeshData::calcCellCells
    (
        mesh,
        agglom,
        agglomPoints.size(),
        true,        // Global mesh connectivity
        cellCells
    );

    labelList decomp;
    decomposeGeneral
    (
        cellCells.values(),
        cellCells.offsets(),
        agglomWeights,
        decomp
    );

    // From coarse back to fine for original mesh
    return labelList(decomp, agglom);
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
) const
{
    if (!cellCentres.empty() && (cellCentres.size() != globalCellCells.size()))
    {
        FatalErrorInFunction
            << "Number of cell centres (" << cellCentres.size()
            << ") != number of cells (" << globalCellCells.size() << ")"
            << exit(FatalError);
    }

    // CompactListList is already
    // Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    labelList decomp;
    decomposeGeneral
    (
        globalCellCells.values(),
        globalCellCells.offsets(),
        cellWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
) const
{
    if (!cellCentres.empty() && (cellCentres.size() != globalCellCells.size()))
    {
        FatalErrorInFunction
            << "Number of cell centres (" << cellCentres.size()
            << ") != number of cells (" << globalCellCells.size() << ")"
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli


    auto cellCells(CompactListList<label>::pack(globalCellCells));

    labelList decomp;
    decomposeGeneral
    (
        cellCells.values(),
        cellCells.offsets(),
        cellWeights,
        decomp
    );

    return decomp;
}

// ************************************************************************* //
