/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "multiLevelDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiLevelDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiLevelDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiLevelDecomp::createMethodsDict()
{
    methodsDict_.clear();

    word defaultMethod;
    labelList domains;

    label nTotal = 0;
    label nLevels = 0;

    // Found (non-recursive, no patterns) "method" and "domains" ?
    // Allow as quick short-cut entry
    if
    (
        // non-recursive, no patterns
        coeffsDict_.readIfPresent("method", defaultMethod, keyType::LITERAL)
        // non-recursive, no patterns
     && coeffsDict_.readIfPresent("domains", domains, keyType::LITERAL)
    )
    {
        // Short-cut version specified by method, domains only

        nTotal = (domains.empty() ? 0 : 1);

        for (const label n : domains)
        {
            nTotal *= n;
            ++nLevels;
        }

        if (nTotal == 1)
        {
            // Emit Warning
            nTotal = nDomains();
            nLevels = 1;

            domains.setSize(1);
            domains[0] = nTotal;
        }
        else if (nTotal > 0 && nTotal < nDomains() && !(nDomains() % nTotal))
        {
            // nTotal < nDomains, but with an integral factor,
            // which we insert as level 0
            ++nLevels;

            labelList old(std::move(domains));

            domains.setSize(old.size()+1);

            domains[0] = nDomains() / nTotal;
            forAll(old, i)
            {
                domains[i+1] = old[i];
            }
            nTotal *= domains[0];

            Info<<"    inferred level0 with " << domains[0]
                << " domains" << nl << nl;
        }

        if (!nLevels || nTotal != nDomains())
        {
            FatalErrorInFunction
                << "Top level decomposition specifies " << nDomains()
                << " domains which is not equal to the product of"
                << " all sub domains " << nTotal
                << exit(FatalError);
        }

        // Create editable methods dictionaries
        nLevels = 0;

        // Common coeffs dictionary
        const dictionary& subMethodCoeffsDict
        (
            findCoeffsDict
            (
                coeffsDict_,
                defaultMethod + "Coeffs",
                selectionType::NULL_DICT
            )
        );

        for (const label n : domains)
        {
            const word levelName("level" + Foam::name(nLevels++));

            entry* dictptr = methodsDict_.set(levelName, dictionary());

            dictionary& dict = dictptr->dict();
            dict.add("method", defaultMethod);
            dict.add("numberOfSubdomains", n);

            // Inject coeffs dictionary too
            if (subMethodCoeffsDict.size())
            {
                dict.add(subMethodCoeffsDict.dictName(), subMethodCoeffsDict);
            }
        }
    }
    else
    {
        // Specified by full dictionaries

        // Create editable methods dictionaries
        // - Only consider sub-dictionaries with a "numberOfSubdomains" entry
        //   This automatically filters out any coeffs dictionaries

        for (const entry& dEntry : coeffsDict_)
        {
            word methodName;

            if
            (
                dEntry.isDict()
                // non-recursive, no patterns
             && dEntry.dict().found("numberOfSubdomains", keyType::LITERAL)
            )
            {
                // No method specified? can use a default method?

                const bool addDefaultMethod
                (
                    !(dEntry.dict().found("method", keyType::LITERAL))
                 && !defaultMethod.empty()
                );

                entry* e = methodsDict_.add(dEntry);

                if (addDefaultMethod && e && e->isDict())
                {
                    e->dict().add("method", defaultMethod);
                }
            }
        }
    }
}


void Foam::multiLevelDecomp::setMethods()
{
    // Assuming methodsDict_ has be properly created, convert the method
    // dictionaries to actual methods

    label nLevels = 0;

    methods_.resize_null(methodsDict_.size());
    for (const entry& dEntry : methodsDict_)
    {
        // Dictionary entries only
        // - these method dictionaries are non-regional
        if (dEntry.isDict())
        {
            methods_.set
            (
                nLevels++,
                // non-verbose would be nicer
                decompositionMethod::New(dEntry.dict())
            );
        }
    }

    methods_.resize(nLevels);

    // Verify that nTotal is correct based on what each method delivers

    Info<< nl
        << "Decompose " << type() << " [" << nDomains() << "] in "
        << nLevels << " levels:" << endl;

    label nTotal = 1;
    forAll(methods_, i)
    {
        Info<< "    level " << i << " : " << methods_[i].type()
            << " [" << methods_[i].nDomains() << "]" << endl;

        nTotal *= methods_[i].nDomains();
    }

    if (nTotal != nDomains())
    {
        FatalErrorInFunction
            << "Top level decomposition specifies " << nDomains()
            << " domains which is not equal to the product of"
            << " all sub domains " << nTotal
            << exit(FatalError);
    }
}


// Given a subset of cells determine the new global indices. The problem
// is in the cells from neighbouring processors which need to be renumbered.
void Foam::multiLevelDecomp::subsetGlobalCellCells
(
    const label nDomains,
    const label domainI,
    const labelList& dist,

    const labelListList& cellCells,
    const labelList& set,
    labelListList& subCellCells,
    labelList& cutConnections
) const
{
    const globalIndex globalCells(cellCells.size());

    // Determine new index for cells by inverting subset
    labelList oldToNew(invert(cellCells.size(), set));

    // Subset locally the elements for which I have data
    subCellCells = UIndirectList<labelList>(cellCells, set);

    // Get new indices for neighbouring processors
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, subCellCells, compactMap);
    map.distribute(oldToNew);
    labelList allDist(dist);
    map.distribute(allDist);

    // Now we have:
    // oldToNew : the locally-compact numbering of all our cellCells. -1 if
    //            cellCell is not in set.
    // allDist  : destination domain for all our cellCells
    // subCellCells : indexes into oldToNew and allDist

    // Globally compact numbering for cells in set.
    const globalIndex globalSubCells(set.size());

    // Now subCellCells contains indices into oldToNew which are the
    // new locations of the neighbouring cells.

    cutConnections.resize_nocopy(nDomains);
    cutConnections = 0;

    forAll(subCellCells, subCelli)
    {
        labelList& cCells = subCellCells[subCelli];

        // Keep the connections to valid mapped cells
        label newI = 0;
        forAll(cCells, i)
        {
            // Get locally-compact cell index of neighbouring cell
            const label nbrCelli = oldToNew[cCells[i]];
            if (nbrCelli == -1)
            {
                cutConnections[allDist[cCells[i]]]++;
            }
            else
            {
                // Reconvert local cell index into global one

                // Get original neighbour
                const label celli = set[subCelli];
                const label oldNbrCelli = cellCells[celli][i];
                // Get processor from original neighbour
                const label proci = globalCells.whichProcID(oldNbrCelli);
                // Convert into global compact numbering
                cCells[newI++] = globalSubCells.toGlobal(proci, nbrCelli);
            }
        }
        cCells.setSize(newI);
    }
}


void Foam::multiLevelDecomp::decompose
(
    const labelListList& pointPoints,
    const pointField& points,
    const scalarField& pointWeights,
    const labelUList& pointMap,     // map back to original points
    const label currLevel,
    const label leafOffset,

    labelList& finalDecomp
) const
{
    labelList dist
    (
        methods_[currLevel].decompose
        (
            pointPoints,
            points,
            pointWeights
        )
    );

    // The next recursion level
    const label nextLevel = currLevel+1;

    // Number of domains at this current level
    const label nCurrDomains = methods_[currLevel].nDomains();

    // Calculate the domain remapping.
    // The decompose() method delivers a distribution of [0..nDomains-1]
    // which we map to the final location according to the decomposition
    // leaf we are on.

    labelList domainLookup(nCurrDomains);
    {
        label sizes = 1;  // Cumulative number of domains
        for (label i = 0; i <= currLevel; ++i)
        {
            sizes *= methods_[i].nDomains();
        }

        // Distribution of domains at this level
        sizes = this->nDomains() / sizes;

        forAll(domainLookup, i)
        {
            domainLookup[i] = i * sizes + leafOffset;
        }
    }

    if (debug)
    {
        Info<< "Distribute at level " << currLevel
            << " to domains" << nl
            << flatOutput(domainLookup) << endl;
    }

    // Extract processor+local index from point-point addressing
    forAll(pointMap, i)
    {
        const label orig = pointMap[i];
        finalDecomp[orig] = domainLookup[dist[i]];
    }

    if (nextLevel < methods_.size())
    {
        // Recurse

        // Extract processor+local index from point-point addressing
        if (debug && Pstream::master())
        {
            Pout<< "Decomposition at level " << currLevel << " :" << endl;
        }

        for (label domainI = 0; domainI < nCurrDomains; ++domainI)
        {
            // Extract elements for current domain
            const labelList domainPoints(findIndices(dist, domainI));

            // Subset point-wise data.
            pointField subPoints(points, domainPoints);
            scalarField subWeights;
            if (pointWeights.size() == points.size())
            {
                subWeights = scalarField(pointWeights, domainPoints);
            }

            labelList subPointMap(labelUIndList(pointMap, domainPoints));
            // Subset point-point addressing (adapt global numbering)
            labelListList subPointPoints;
            labelList nOutsideConnections;
            subsetGlobalCellCells
            (
                nCurrDomains,
                domainI,
                dist,

                pointPoints,
                domainPoints,

                subPointPoints,
                nOutsideConnections
            );

            label nPoints = returnReduce(domainPoints.size(), sumOp<label>());

            Pstream::listCombineReduce(nOutsideConnections, plusEqOp<label>());
            label nPatches = 0;
            label nFaces = 0;
            for (const label nConnect : nOutsideConnections)
            {
                if (nConnect > 0)
                {
                    ++nPatches;
                    nFaces += nConnect;
                }
            }

            string oldPrefix;
            if (debug && Pstream::master())
            {
                Pout<< "    Domain " << domainI << nl
                    << "        Number of cells = " << nPoints << nl
                    << "        Number of inter-domain patches = " << nPatches
                    << nl
                    << "        Number of inter-domain faces = " << nFaces << nl
                    << endl;
                oldPrefix = Pout.prefix();
                Pout.prefix() = "  " + oldPrefix;
            }

            decompose
            (
                subPointPoints,
                subPoints,
                subWeights,
                subPointMap,
                nextLevel,
                domainLookup[domainI], // The offset for this level and leaf

                finalDecomp
            );
            if (debug && Pstream::master())
            {
                Pout.prefix() = oldPrefix;
            }
        }


        if (debug)
        {
            // Do straight decompose of two levels
            const label nNext = methods_[nextLevel].nDomains();
            const label nTotal = nCurrDomains * nNext;

            // Get original level0 dictionary and modify numberOfSubdomains
            dictionary level0Dict;
            for (const entry& dEntry : methodsDict_)
            {
                if (dEntry.isDict())
                {
                    level0Dict = dEntry.dict();
                    break;
                }
            }
            level0Dict.set("numberOfSubdomains", nTotal);

            if (debug && Pstream::master())
            {
                Pout<< "Reference decomposition with " << level0Dict << " :"
                    << endl;
            }

            autoPtr<decompositionMethod> method0 = decompositionMethod::New
            (
                level0Dict
            );
            labelList dist
            (
                method0().decompose
                (
                    pointPoints,
                    points,
                    pointWeights
                )
            );

            for (label blockI = 0; blockI < nCurrDomains; ++blockI)
            {
                // Count the number inbetween blocks of nNext size

                label nPoints = 0;
                labelList nOutsideConnections(nCurrDomains, Foam::zero{});
                forAll(pointPoints, pointi)
                {
                    if ((dist[pointi] / nNext) == blockI)
                    {
                        nPoints++;

                        const labelList& pPoints = pointPoints[pointi];

                        forAll(pPoints, i)
                        {
                            const label distBlockI = dist[pPoints[i]] / nNext;
                            if (distBlockI != blockI)
                            {
                                nOutsideConnections[distBlockI]++;
                            }
                        }
                    }
                }

                reduce(nPoints, sumOp<label>());
                Pstream::listCombineReduce
                (
                    nOutsideConnections,
                    plusEqOp<label>()
                );

                label nPatches = 0;
                label nFaces = 0;
                for (const label nConnect : nOutsideConnections)
                {
                    if (nConnect > 0)
                    {
                        ++nPatches;
                        nFaces += nConnect;
                    }
                }

                if (debug && Pstream::master())
                {
                    Pout<< "    Domain " << blockI << nl
                        << "        Number of cells = " << nPoints << nl
                        << "        Number of inter-domain patches = "
                        << nPatches << nl
                        << "        Number of inter-domain faces = " << nFaces
                        << nl << endl;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiLevelDecomp::multiLevelDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_
    (
        findCoeffsDict
        (
            typeName + "Coeffs",
            (selectionType::EXACT | selectionType::MANDATORY)
        )
    ),
    methodsDict_(),
    methods_()
{
    createMethodsDict();
    setMethods();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiLevelDecomp::parallelAware() const
{
    for (const decompositionMethod& meth : methods_)
    {
        if (!meth.parallelAware())
        {
            return false;
        }
    }

    return true;
}


Foam::labelList Foam::multiLevelDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& cc,
    const scalarField& cWeights
) const
{
    labelList finalDecomp(cc.size(), Foam::zero{});
    labelList cellMap(Foam::identity(cc.size()));

    CompactListList<label> cellCells;
    globalMeshData::calcCellCells
    (
        mesh,
        cellMap,
        cellMap.size(),
        true,       // Global mesh connectivity
        cellCells
    );

    decompose
    (
        cellCells.unpack(),
        cc,
        cWeights,
        cellMap,      // map back to original cell points
        0,            // currLevel
        0,            // leafOffset

        finalDecomp
    );

    return finalDecomp;
}


Foam::labelList Foam::multiLevelDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights
) const
{
    labelList finalDecomp(cc.size(), Foam::zero{});
    labelList cellMap(Foam::identity(cc.size()));

    decompose
    (
        globalCellCells.unpack(),
        cc,
        cWeights,
        cellMap,      // map back to original cell points
        0,            // currLevel
        0,            // leafOffset

        finalDecomp
    );

    return finalDecomp;
}


Foam::labelList Foam::multiLevelDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights
) const
{
    labelList finalDecomp(cc.size(), Foam::zero{});
    labelList cellMap(Foam::identity(cc.size()));

    decompose
    (
        globalCellCells,
        cc,
        cWeights,
        cellMap,      // map back to original cell points
        0,            // currLevel
        0,            // leafOffset

        finalDecomp
    );

    return finalDecomp;
}


// ************************************************************************* //
