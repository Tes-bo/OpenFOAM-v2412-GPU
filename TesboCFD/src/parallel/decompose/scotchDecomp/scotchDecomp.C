/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "Time.H"
#include "PrecisionAdaptor.H"
#include "OFstream.H"
#include <cstdio>
#include <limits>

// Probably not needed, but in case we pickup a ptscotch.h ...
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include "scotch.h"

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
    defineTypeNameAndDebug(scotchDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::scotchDecomp::decomposeSerial
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


    // Dump graph
    if (coeffsDict_.getOrDefault("writeGraph", false))
    {
        OFstream str(graphPath_ + ".grf");

        Info<< "Dumping Scotch graph file to " << str.name() << nl
            << "Use this in combination with gpart." << endl;

        const label numConnect = adjncy.size();

        // Version 0 = Graph file (.grf)
        str << "0" << nl;

        // Number of vertices,
        // number of edges (connections)
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
        DebugInfo << "scotchDecomp : Using strategy " << strategy << endl;

        SCOTCH_stratGraphMap(&stradat, strategy.c_str());
        //fprintf(stdout, "S\tStrat=");
        //SCOTCH_stratSave(&stradat, stdout);
        //fprintf(stdout, "\n");
    }


    // Graph
    // ~~~~~

    // Check for externally provided cellweights and if so initialise weights

    bool hasWeights = !cWeights.empty();

    // Note: min, not gMin since routine runs on master only.
    const scalar minWeights = hasWeights ? min(cWeights) : scalar(1);

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

        const scalar velotabSum = sum(cWeights)/minWeights;

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
    }
    else
    {
        // HACK: specify uniform weights
        // - seems that scotch takes different code paths internally
        //   if weights are not specified (issue #3063)

        velotab.resize(numCells);
        velotab = static_cast<SCOTCH_Num>(1);
    }


    //
    // Decomposition graph
    //

    SCOTCH_Graph grafdat;
    check
    (
        SCOTCH_graphInit(&grafdat),
        "SCOTCH_graphInit"
    );
    check
    (
        SCOTCH_graphBuild
        (
            &grafdat,               // Graph to build
            0,                      // Base for indexing (C-style)

            numCells,               // Number of vertices [== nCells]
            xadj_param().cdata(),   // verttab, start index per cell into adjncy
            nullptr,                // vendtab, end index (nullptr == automatic)

            velotab.cdata(),        // velotab, vertex weights
            nullptr,                // Vertex labels (nullptr == ignore)

            adjncy.size(),          // Number of graph edges
            adjncy_param().cdata(), // Edge array
            nullptr                 // Edge weights (nullptr == ignore)
        ),
        "SCOTCH_graphBuild"
    );
    check
    (
        SCOTCH_graphCheck(&grafdat),
        "SCOTCH_graphCheck"
    );


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check
    (
        SCOTCH_archInit(&archdat),
        "SCOTCH_archInit"
    );

    List<SCOTCH_Num> procWeights;

    // Do optional tree-like/multi-level decomposition by specifying
    // domains/weights
    List<SCOTCH_Num> domains;
    List<scalar> dWeights;

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
            << "scotchDecomp : Using procesor weights "
            << procWeights << endl;

        check
        (
            SCOTCH_archCmpltw(&archdat, nDomains_, procWeights.cdata()),
            "SCOTCH_archCmpltw"
        );
    }
    else if
    (
        coeffsDict_.readIfPresent("domains", domains, keyType::LITERAL)
     && coeffsDict_.readIfPresent("domainWeights", dWeights, keyType::LITERAL)
    )
    {
        // multi-level

        label nTotal = 1;
        for (const label n : domains)
        {
            nTotal *= n;
        }

        if (nTotal < nDomains())
        {
            const label sz = domains.size();
            domains.setSize(sz+1);
            dWeights.setSize(sz+1);
            for (label i = sz-1; i >= 0; i--)
            {
                domains[i+1] = domains[i];
                dWeights[i+1] = dWeights[i];
            }

            if (nDomains() % nTotal)
            {
                FatalErrorInFunction
                    << "Top level decomposition specifies " << nDomains()
                    << " domains which is not equal to the product of"
                    << " all sub domains " << nTotal
                    << exit(FatalError);
            }

            domains[0] = nDomains() / nTotal;
            dWeights[0] = scalar(1);
        }


        // Note: min, not gMin since routine runs on master only.
        const scalar minWeights = min(dWeights);

        // Convert to integers.
        List<SCOTCH_Num> weights(dWeights.size());

        forAll(weights, i)
        {
            weights[i] = static_cast<SCOTCH_Num>
            (
                (dWeights[i]/minWeights - 1) + 1
            );
        }

        check
        (
            SCOTCH_archTleaf
            (
                &archdat,
                SCOTCH_Num(domains.size()),
                domains.cdata(),
                weights.cdata()
            ),
            "SCOTCH_archTleaf"
        );
    }
    else
    {
        check
        (
            SCOTCH_archCmplt(&archdat, nDomains_),
            "SCOTCH_archCmplt"
        );


        //- Hack to test clustering. Note that decomp is non-compact
        //  numbers!
        //
        ////- Set up variable sizes architecture
        //check
        //(
        //    SCOTCH_archVcmplt(&archdat),
        //    "SCOTCH_archVcmplt"
        //);
        //
        ////- Stategy flags: go for quality or load balance (or leave default)
        //SCOTCH_Num straval = 0;
        ////straval |= SCOTCH_STRATQUALITY;
        ////straval |= SCOTCH_STRATQUALITY;
        //
        ////- Number of cells per agglomeration
        ////SCOTCH_Num agglomSize = SCOTCH_archSize(&archdat);
        //SCOTCH_Num agglomSize = 3;
        //
        ////- Build strategy for agglomeration
        //check
        //(
        //    SCOTCH_stratGraphClusterBuild
        //    (
        //        &stradat,   // strategy to build
        //        straval,    // strategy flags
        //        agglomSize, // cells per cluster
        //        1.0,        // weight?
        //        0.01        // max load imbalance
        //    ),
        //    "SCOTCH_stratGraphClusterBuild"
        //);
    }


    //SCOTCH_Mapping mapdat;
    //SCOTCH_graphMapInit(&grafdat, &mapdat, &archdat, nullptr);
    //SCOTCH_graphMapCompute(&grafdat, &mapdat, &stradat); /* Perform mapping */
    //SCOTCH_graphMapExit(&grafdat, &mapdat);


    // Hack:switch off fpu error trapping
    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    check
    (
        SCOTCH_graphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            decomp_param.ref().data() // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    //decomp.resize(numCells);
    //check
    //(
    //    SCOTCH_graphPart
    //    (
    //        &grafdat,
    //        nDomains_,           // partnbr
    //        &stradat,            // const SCOTCH_Strat *
    //        decomp_param.ref().data()  // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    SCOTCH_graphExit(&grafdat);     // Release storage for graph
    SCOTCH_stratExit(&stradat);     // Release storage for strategy
    SCOTCH_archExit(&archdat);      // Release storage for network topology

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp(const label numDomains)
:
    metisLikeDecomp(numDomains)
{}


Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    // Where to write graph
    graphPath_ = getGraphPathBase(mesh);

    return metisLikeDecomp::decompose
    (
        mesh,
        points,
        pointWeights
    );
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
) const
{
    // Where to write graph
    graphPath_ = getGraphPathBase(mesh);

    return metisLikeDecomp::decompose
    (
        mesh,
        agglom,
        agglomPoints,
        agglomWeights
    );
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    // Where to write graph
    graphPath_ = "scotch.grf";

    return metisLikeDecomp::decompose
    (
        globalCellCells,
        cellCentres,
        cWeights
    );
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    // Where to write graph
    graphPath_ = "scotch.grf";

    return metisLikeDecomp::decompose
    (
        globalCellCells,
        cellCentres,
        cWeights
    );
}


// ************************************************************************* //
