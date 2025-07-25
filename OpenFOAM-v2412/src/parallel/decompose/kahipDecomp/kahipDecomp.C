/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "kahipDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "PrecisionAdaptor.H"
#include <map>
#include <string>
#include <vector>

// Probably not needed, but in case we pickup parhip_interface.h
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include "kaHIP_interface.h"

// Provide a clear error message if we have a severe size mismatch
// Allow widening, but not narrowing

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kahipDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        kahipDecomp,
        dictionary
    );
}


const Foam::Enum
<
    Foam::kahipDecomp::configs
>
Foam::kahipDecomp::configNames
({
    { kahipDecomp::configs::FAST, "fast" },
    { kahipDecomp::configs::ECO, "eco" },
    { kahipDecomp::configs::STRONG, "strong" },
    { kahipDecomp::configs::FASTSOCIAL, "fast-social" },
    { kahipDecomp::configs::ECOSOCIAL, "eco-social" },
    { kahipDecomp::configs::STRONGSOCIAL, "strong-social" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::kahipDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    // Default setup
    enum configs kahipConfig = configs::FAST;
    double imbalance = 0.01;
    int seed = 0;
    bool verbose = false;

    #if WM_LABEL_SIZE == 64
    if (xadj.size()-1 > INT_MAX)
    {
        FatalErrorInFunction
            << "Cannot decompose " << (xadj.size()-1) << " cells," << nl
            << "Exceeded integer limit of " << INT_MAX << nl
            << exit(FatalError);
    }
    #endif

    int numCells = Foam::max(0, (xadj.size()-1));

    // Addressing
    ConstPrecisionAdaptor<int, label, List> adjncy_param(adjncy);
    ConstPrecisionAdaptor<int, label, List> xadj_param(xadj);

    // Output: cell -> processor addressing
    decomp.resize_nocopy(numCells);
    decomp = 0;
    PrecisionAdaptor<int, label, List> decomp_param(decomp, false);

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
            << "Number of weights (" << cWeights.size()
            << ") != number of cells (" << numCells << ")"
            << exit(FatalError);
    }

    // Cell weights (so on the vertices of the dual)
    List<int> cellWeights;

    if (hasWeights)
    {
        // Convert to integers.
        cellWeights.resize_nocopy(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = static_cast<int>
            (
                cWeights[i]/minWeights
            );
        }
    }

    configNames.readIfPresent("config", coeffsDict_, kahipConfig);
    coeffsDict_.readIfPresent("imbalance", imbalance);
    coeffsDict_.readIfPresent("verbose", verbose);

    Info<< "kahipDecomp :"
        << " config=" << configNames[kahipConfig]
        << " imbalance=" << imbalance;

    if (coeffsDict_.readIfPresent("seed", seed))
    {
        Info<< " seed=" << seed;
    }

    // Additional sizing parameters (testing only)
    std::map<std::string, std::vector<int>> sizingParams;

    List<int> labels;
    if
    (
        coeffsDict_.readIfPresent("hierarchy", labels)
     && !labels.empty()
    )
    {
        std::vector<int> vec;
        vec.reserve(labels.size()+1);

        // Verify sizing

        int n = 1;
        for (const auto val : labels)
        {
            n *= val;
            vec.push_back(val);
        }

        if (n != nDomains_)
        {
            // Size mismatch. Try to correct.

            if (nDomains_ % n)
            {
                WarningInFunction
                    << "Mismatch in number of processors and "
                    << "hierarchy specified" << flatOutput(labels) << endl;

                vec.clear();
            }
            else
            {
                // Evenly divisible, add extra hierarchy level
                vec.push_back(nDomains_ / n);
            }
        }

        if (!vec.empty())
        {
            sizingParams["hierarchy"] = std::move(vec);
            Info<< " hierarchy=" << flatOutput(labels);
        }
    }

    if
    (
        coeffsDict_.readIfPresent("distance", labels)
     && !labels.empty()
    )
    {
        std::vector<int> vec(labels.size());

        forAll(labels, i)
        {
            vec[i] = labels[i];
        }

        sizingParams["distance"] = std::move(vec);
        Info<< " distance=" << flatOutput(labels);
    }

    Info<< endl;


    // Number of partitions
    int nParts = nDomains_;

    // Output: number of cut edges
    int edgeCut = 0;


#if 0 // WIP: #ifdef KAFFPA_CPP_INTERFACE
    kaffpa_cpp
    (
        &numCells,          // num vertices in graph
        (cellWeights.empty() ? nullptr : cellWeights.data()), // vertex wts
        xadj_param.constCast().data(),          // indexing into adjncy
        nullptr,            // edge wts
        adjncy_param.constCast().data(),        // neighbour info
        &nParts,            // nparts
        &imbalance,         // amount of imbalance allowed
        !verbose,           // suppress output
        seed,               // for random
        int(kahipConfig),
        &edgeCut,                   // [output]
        decomp_param.ref().data(),  // [output]
        sizingParams
    );
#else
    kaffpa
    (
        &numCells,          // num vertices in graph
        (cellWeights.empty() ? nullptr : cellWeights.data()), // vertex wts
        xadj_param.constCast().data(),          // indexing into adjncy
        nullptr,            // edge wts
        adjncy_param.constCast().data(),        // neighbour info
        &nParts,            // nparts
        &imbalance,         // amount of imbalance allowed
        !verbose,           // suppress output
        seed,               // for random
        int(kahipConfig),
        &edgeCut,                   // [output]
        decomp_param.ref().data()   // [output]
    );
#endif

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kahipDecomp::kahipDecomp(const label numDomains)
:
    metisLikeDecomp(numDomains)
{}


Foam::kahipDecomp::kahipDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// ************************************************************************* //
