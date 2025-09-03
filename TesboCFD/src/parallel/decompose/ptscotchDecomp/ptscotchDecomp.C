/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "ptscotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "globalMeshData.H"
#include "Time.H"
#include "PrecisionAdaptor.H"
#include "OFstream.H"
#include <cstdio>
#include <limits>

// Include MPI without any C++ bindings
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

// Avoid too many warnings from mpi.h
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <mpi.h>
#include "ptscotch.h"

// Hack: scotch generates floating point errors so need to switch off error
//       trapping!
#ifdef __GLIBC__
    #ifndef _GNU_SOURCE
        #define _GNU_SOURCE
    #endif
    #include <fenv.h>
#endif

// Error if we attempt narrowing
static_assert
(
    sizeof(Foam::label) <= sizeof(SCOTCH_Num),
    "SCOTCH_Num is too small for Foam::label, check your scotch headers"
);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ptscotchDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotchDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check and print error message
static inline void check(const int retVal, const char* what)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << what
            << " failed (" << retVal << ")\n"
            << exit(FatalError);
    }
}

// The mesh-relative graph path/name (without extension)
static inline Foam::fileName getGraphPathBase(const polyMesh& mesh)
{
    return mesh.time().path()/mesh.name();
}


} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::ptscotchDecomp::decompose
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    const SCOTCH_Num numCells = Foam::max(0, (xadj.size()-1));

    // Addressing
    ConstPrecisionAdaptor<SCOTCH_Num, label, List> adjncy_param(adjncy);
    ConstPrecisionAdaptor<SCOTCH_Num, label, List> xadj_param(xadj);

    // Output: cell -> processor addressing
    decomp.resize_nocopy(numCells);
    decomp = 0;
    PrecisionAdaptor<SCOTCH_Num, label, List> decomp_param(decomp, false);

    // Avoid potential nullptr issues with zero-sized arrays
    labelList adjncy_dummy, xadj_dummy, decomp_dummy;
    if (!numCells)
    {
        adjncy_dummy.resize(1, 0);
        adjncy_param.set(adjncy_dummy);

        xadj_dummy.resize(2, 0);
        xadj_param.set(xadj_dummy);

        decomp_dummy.resize(1, 0);
        decomp_param.clear();  // Avoid propagating spurious values
        decomp_param.set(decomp_dummy);
    }


    if (debug & 2)
    {
        Pout<< "ptscotchDecomp : " << numCells << " cells" << endl;
    }

    // Dump graph
    if (coeffsDict_.getOrDefault("writeGraph", false))
    {
        OFstream str
        (
            graphPath_ + "_" + Foam::name(Pstream::myProcNo()) + ".dgr"
        );

        Pout<< "Dumping Scotch graph file to " << str.name() << endl
            << "Use this in combination with dgpart." << endl;

        const label numConnect = adjncy.size();
        const label nTotCells = returnReduce(numCells, sumOp<label>());
        const label nTotConnect = returnReduce(numConnect, sumOp<label>());

        // Version 2 = Distributed graph file (.dgrf)
        str << "2" << nl;

        // Number of files (procglbnbr), my file number (procloc)
        str << Pstream::nProcs() << ' ' << Pstream::myProcNo() << nl;

        // Total number of vertices (vertglbnbr),
        // Total number of connections (edgeglbnbr)
        str << nTotCells << ' ' << nTotConnect << nl;

        // Local number of vertices (vertlocnbr),
        // Local number of connections (edgelocnbr)
        str << numCells << ' ' << numConnect << nl;

        // Numbering starts from 0
        // 100*hasVertlabels+10*hasEdgeWeights+1*hasVertWeights
        str << "0 000" << nl;

        for (label celli = 0; celli < numCells; ++celli)
        {
            const label beg = xadj[celli];
            const label end = xadj[celli+1];

            str << (end-beg);  // size

            for (label i = beg; i < end; ++i)
            {
                str << ' ' << adjncy[i];
            }
            str << nl;
        }
    }


    // Make repeatable
    SCOTCH_randomReset();

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check
    (
        SCOTCH_stratInit(&stradat),
        "SCOTCH_stratInit"
    );

    string strategy;
    if (coeffsDict_.readIfPresent("strategy", strategy))
    {
        DebugInfo
            << "ptscotchDecomp : Using strategy " << strategy << endl;

        SCOTCH_stratDgraphMap(&stradat, strategy.c_str());
        //fprintf(stdout, "S\tStrat=");
        //SCOTCH_stratSave(&stradat, stdout);
        //fprintf(stdout, "\n");
    }


    // Graph
    // ~~~~~

    // Check for externally provided cellweights and if so initialise weights

    bool hasWeights = returnReduceOr(cWeights.size());

    const scalar minWeights = hasWeights ? gMin(cWeights) : scalar(1);

    if (minWeights <= 0)
    {
        hasWeights = false;
        WarningInFunction
            << "Illegal minimum weight " << minWeights
            << " ... ignoring"
            << endl;
    }
    else if (hasWeights && (cWeights.size() != numCells))
    {
        FatalErrorInFunction
            << "Number of cell weights " << cWeights.size()
            << " does not equal number of cells " << numCells
            << exit(FatalError);
    }


    List<SCOTCH_Num> velotab;

    if (hasWeights)
    {
        scalar rangeScale(1);

        const scalar velotabSum = gSum(cWeights)/minWeights;

        const scalar upperRange = static_cast<scalar>
        (
            std::numeric_limits<SCOTCH_Num>::max()-1
        );

        if (velotabSum > upperRange)
        {
            // 0.9 factor of safety to avoid floating point round-off in
            // rangeScale tipping the subsequent sum over the integer limit.
            rangeScale = 0.9*upperRange/velotabSum;

            WarningInFunction
                << "Sum of weights overflows SCOTCH_Num: " << velotabSum
                << ", compressing by factor " << rangeScale << endl;
        }

        if (cWeights.size())
        {
            // Convert to integers.
            velotab.resize(cWeights.size());

            forAll(velotab, i)
            {
                velotab[i] = static_cast<SCOTCH_Num>
                (
                    ((cWeights[i]/minWeights - 1)*rangeScale) + 1
                );
            }
        }
        else
        {
            // Locally zero cells but not globally.
            // Provide some size to avoid null pointer.
            velotab.resize(1, 1);
        }
    }
    else
    {
        // HACK: specify uniform weights
        // - seems that scotch takes different code paths internally
        //   if weights are not specified (issue #3063)

        if (numCells > 0)
        {
            velotab.resize(numCells);
        }
        else
        {
            // Locally zero cells but not globally.
            // Provide some size to avoid null pointer.
            velotab.resize(1);
        }

        velotab = static_cast<SCOTCH_Num>(1);
    }


    //
    // Decomposition graph
    //

    if (debug & 2)
    {
        Pout<< "SCOTCH_dgraphInit" << endl;
    }
    SCOTCH_Dgraph grafdat;
    check
    (
        SCOTCH_dgraphInit(&grafdat, MPI_COMM_WORLD),
        "SCOTCH_dgraphInit"
    );

    if (debug & 2)
    {
        Pout<< "SCOTCH_dgraphBuild with:" << nl
            << "numCells        : " << numCells << nl
            << "xadj            : " << name(xadj_param().cdata()) << nl
            << "velotab         : " << name(velotab.cdata()) << nl
            << "adjncySize      : " << adjncy_param().size() << nl
            << "adjncy          : " << name(adjncy_param().cdata()) << nl
            << endl;
    }

    check
    (
        SCOTCH_dgraphBuild
        (
            &grafdat,               // Graph to build
            0,                      // Base for indexing (C-style)

            numCells,               // vertlocnbr [== nCells]
            numCells,               // vertlocmax

            xadj_param.constCast().data(),
                                    // vertloctab, start index per cell into
                                    // adjncy
            (xadj_param.constCast().data()+1),
                                    // vendloctab, end index  ,,

            velotab.data(),         // veloloctab, vtx weights
            nullptr,                // vlblloctab

            adjncy.size(),          // edgelocnbr, number of arcs
            adjncy.size(),          // edgelocsiz
            adjncy_param.constCast().data(), // edgeloctab
            nullptr,                // edgegsttab
            nullptr                 // edlotab, edge weights
        ),
        "SCOTCH_dgraphBuild"
    );

    if (debug & 2)
    {
        Pout<< "SCOTCH_dgraphCheck" << endl;
    }
    check
    (
        SCOTCH_dgraphCheck(&grafdat),
        "SCOTCH_dgraphCheck"
    );


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    if (debug & 2)
    {
        Pout<< "SCOTCH_archInit" << endl;
    }
    SCOTCH_Arch archdat;
    check
    (
        SCOTCH_archInit(&archdat),
        "SCOTCH_archInit"
    );

    List<SCOTCH_Num> procWeights;
    if
    (
        coeffsDict_.readIfPresent("processorWeights", procWeights)
     && !procWeights.empty()
    )
    {
        if (procWeights.size() != nDomains_)
        {
            FatalIOErrorInFunction(coeffsDict_)
                << "processorWeights (" << procWeights.size()
                << ") != number of domains (" << nDomains_ << ")" << nl
                << exit(FatalIOError);
        }

        DebugInfo
            << "ptscotchDecomp : Using procesor weights "
            << procWeights << endl;

        check
        (
            SCOTCH_archCmpltw(&archdat, nDomains_, procWeights.cdata()),
            "SCOTCH_archCmpltw"
        );
    }
    else
    {
        // Check for optional domains/weights
        List<SCOTCH_Num> domains;
        List<scalar> dWeights;
        if
        (
            coeffsDict_.readIfPresent
            (
                "domains",
                domains,
                keyType::LITERAL
            )
         && coeffsDict_.readIfPresent
            (
                "domainWeights",
                dWeights,
                keyType::LITERAL
            )
        )
        {
            WarningInFunction
                << "Ignoring multi-level decomposition since"
                << " not supported by ptscotch."
                << " It is supported by scotch" << endl;
        }

        if (debug & 2)
        {
            Pout<< "SCOTCH_archCmplt" << endl;
        }
        check
        (
            SCOTCH_archCmplt(&archdat, nDomains_),
            "SCOTCH_archCmplt"
        );
    }


    // Hack:switch off fpu error trapping
    #ifdef  FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    if (debug & 2)
    {
        Pout<< "SCOTCH_dgraphMap" << endl;
    }
    check
    (
        SCOTCH_dgraphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            decomp_param.ref().data() // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef  FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    //check
    //(
    //    SCOTCH_dgraphPart
    //    (
    //        &grafdat,
    //        nDomains_,          // partnbr
    //        &stradat,           // const SCOTCH_Strat *
    //        decomp_param.ref().data() // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    if (debug & 2)
    {
        Pout<< "SCOTCH_dgraphExit" << endl;
    }

    SCOTCH_dgraphExit(&grafdat);    // Release storage for graph
    SCOTCH_stratExit(&stradat);     // Release storage for strategy
    SCOTCH_archExit(&archdat);      // Release storage for network topology

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptscotchDecomp::ptscotchDecomp(const label numDomains)
:
    decompositionMethod(numDomains),
    coeffsDict_(dictionary::null)
{}


Foam::ptscotchDecomp::ptscotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_(findCoeffsDict("scotchCoeffs", selectionType::NULL_DICT))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::ptscotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    // Where to write graph
    graphPath_ = getGraphPathBase(mesh);

    if (!points.empty() && (points.size() != mesh.nCells()))
    {
        FatalErrorInFunction
            << "Number of cell centres (" << points.size()
            << ") != number of cells (" << mesh.nCells() << ")"
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    // Global mesh connectivity
    CompactListList<label> cellCells;
    globalMeshData::calcCellCells(mesh, cellCells, true);

    labelList decomp;
    decompose
    (
        cellCells.values(),
        cellCells.offsets(),
        pointWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
) const
{
    // Where to write graph
    graphPath_ = getGraphPathBase(mesh);

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
        true,       // Global mesh connectivity
        cellCells
    );

    // Decompose using weights
    labelList decomp;
    decompose
    (
        cellCells.values(),
        cellCells.offsets(),
        agglomWeights,
        decomp
    );

    // From coarse back to fine for original mesh
    return labelList(decomp, agglom);
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    // Where to write graph
    graphPath_ = "ptscotch";

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
    decompose
    (
        globalCellCells.values(),
        globalCellCells.offsets(),
        cWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    // Where to write graph
    graphPath_ = "ptscotch";

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
    decompose
    (
        cellCells.values(),
        cellCells.offsets(),
        cWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
