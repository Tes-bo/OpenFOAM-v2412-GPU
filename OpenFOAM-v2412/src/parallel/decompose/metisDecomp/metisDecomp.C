/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "PrecisionAdaptor.H"

// Probably not needed...
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include "metis.h"

// Provide a clear error message if we have a severe size mismatch
// Allow widening, but not narrowing
//
// Metis has an 'idx_t' type, but the IDXTYPEWIDTH define is perhaps
// more future-proof?
//#ifdef IDXTYPEWIDTH
//static_assert
//(
//    sizeof(Foam::label) > (IDXTYPEWIDTH/8),
//    "sizeof(Foam::label) > (IDXTYPEWIDTH/8), check your metis headers"
//);
//#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("recursive");

    const dictionary* coeffsDictPtr = decompDict_.findDict("metisCoeffs");

    idx_t numCells = Foam::max(0, (xadj.size()-1));

    // Decomposition options
    List<idx_t> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.data());

    // Processor weights initialised with no size, only used if specified in
    // a file
    Field<real_t> procWeights;

    // Cell weights (so on the vertices of the dual)
    List<idx_t> cellWeights;

    // Face weights (so on the edges of the dual)
    List<idx_t> faceWeights;

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

    if (hasWeights)
    {
        // Convert to integers.
        cellWeights.resize_nocopy(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = idx_t(cWeights[i]/minWeights);
        }
    }


    // Check for user supplied weights and decomp options
    if (coeffsDictPtr)
    {
        const dictionary& coeffDict = *coeffsDictPtr;

        word weightsFile;

        if (coeffDict.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorInFunction
                    << "Method " << method << " in metisCoeffs in dictionary : "
                    << decompDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis method     " << method
                << nl << endl;
        }

        if (coeffDict.readIfPresent("options", options))
        {
            if (options.size() != METIS_NOPTIONS)
            {
                FatalErrorInFunction
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis options     " << options
                << nl << endl;
        }

        if (coeffDict.readIfPresent("processorWeights", procWeights))
        {
            if (procWeights.size() != nDomains_)
            {
                FatalIOErrorInFunction(coeffDict)
                    << "processorWeights (" << procWeights.size()
                    << ") != number of domains (" << nDomains_ << ")" << nl
                    << exit(FatalIOError);
            }

            procWeights /= sum(procWeights);
        }
    }

    idx_t ncon = 1;
    idx_t nProcs = nDomains_;

    // Addressing
    ConstPrecisionAdaptor<idx_t, label, List> xadj_param(xadj);
    ConstPrecisionAdaptor<idx_t, label, List> adjncy_param(adjncy);

    // Output: cell -> processor addressing
    decomp.resize(numCells);
    PrecisionAdaptor<idx_t, label, List> decomp_param(decomp, false);

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


    //
    // Decompose
    //

    // Output: number of cut edges
    idx_t edgeCut = 0;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,                  // num vertices in graph
            &ncon,                      // num balancing constraints
            xadj_param.constCast().data(),      // indexing into adjncy
            adjncy_param.constCast().data(),    // neighbour info
            cellWeights.data(),         // vertex wts
            nullptr,                    // vsize: total communication vol
            faceWeights.data(),         // edge wts
            &nProcs,                    // nParts
            procWeights.data(),         // tpwgts
            nullptr,                    // ubvec: processor imbalance (default)
            options.data(),
            &edgeCut,
            decomp_param.ref().data()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,                  // num vertices in graph
            &ncon,                      // num balancing constraints
            xadj_param.constCast().data(),      // indexing into adjncy
            adjncy_param.constCast().data(),    // neighbour info
            cellWeights.data(),         // vertex wts
            nullptr,                    // vsize: total communication vol
            faceWeights.data(),         // edge wts
            &nProcs,                    // nParts
            procWeights.data(),         // tpwgts
            nullptr,                    // ubvec: processor imbalance (default)
            options.data(),
            &edgeCut,
            decomp_param.ref().data()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp(const label numDomains)
:
    metisLikeDecomp(numDomains)
{}


Foam::metisDecomp::metisDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// ************************************************************************* //
