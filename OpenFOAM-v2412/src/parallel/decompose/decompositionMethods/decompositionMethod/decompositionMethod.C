/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "decompositionMethod.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "syncTools.H"
#include "faceSet.H"
#include "regionSplit.H"
#include "localPointRegion.H"
#include "minData.H"
#include "BitOps.H"
#include "FaceCellWave.H"

// Compatibility (MAY-2014)
#include "preserveBafflesConstraint.H"
#include "preservePatchesConstraint.H"
#include "preserveFaceZonesConstraint.H"
#include "singleProcessorFaceSetsConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionMethod, 0);
    defineRunTimeSelectionTable(decompositionMethod, dictionary);

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Find named coefficents dictionary, or use default "coeffs"
static inline const dictionary* cfindCoeffsDict
(
    const dictionary& dict,
    const word& coeffsName,
    const bool allowDefault
)
{
    const dictionary* dictptr = dict.findDict(coeffsName);
    if (!dictptr && allowDefault)
    {
        dictptr = dict.findDict("coeffs");
    }
    return dictptr;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::decompositionMethod::nDomains
(
    const dictionary& decompDict,
    const word& regionName
)
{
    label nDomainsRegion = 0;
    label nDomainsGlobal = UPstream::nProcs();

    // Allow numberOfSubdomains to be optional in parallel, which allows
    // for missing files on directories that have not yet been created.

    decompDict.readEntry<label>
    (
        "numberOfSubdomains",
        nDomainsGlobal,
        keyType::REGEX,  // keyType::LITERAL?
        (UPstream::parRun() ? IOobject::LAZY_READ : IOobject::MUST_READ)
    );

    if (!regionName.empty())
    {
        const dictionary& regionDict =
            optionalRegionDict(decompDict, regionName);

        if (regionDict.readIfPresent("numberOfSubdomains", nDomainsRegion))
        {
            if (nDomainsRegion >= 1 && nDomainsRegion <= nDomainsGlobal)
            {
                return nDomainsRegion;
            }

            WarningInFunction
                << "Ignoring region [" << regionName
                << "] numberOfSubdomains: " << nDomainsRegion
                << ", using global: " << nDomainsGlobal << nl
                << endl;
        }
    }

    return nDomainsGlobal;
}


const Foam::dictionary& Foam::decompositionMethod::optionalRegionDict
(
    const dictionary& decompDict,
    const word& regionName
)
{
    const dictionary* dictptr = nullptr;
    if
    (
        !regionName.empty()
     && (dictptr = decompDict.findDict("regions")) != nullptr
    )
    {
        dictptr = dictptr->findDict(regionName);
    }
    return (dictptr ? *dictptr : dictionary::null);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::decompositionMethod::constraintCompat(const word& modelType) const
{
    bool usable = decompDict_.found(modelType);
    if (!usable)
    {
        return false;
    }

    for (const auto& item : constraints_)
    {
        if (modelType == item.type())
        {
            usable = false;
            break;
        }
    }

    if (usable)
    {
        Warning
            << nl << "    Using '" << modelType
            << "' constraint specification." << nl;
    }
    else
    {
        Warning
            << nl << "    Ignoring '" << modelType
            << "' constraint specification - was already specified." << nl;
    }

    // The syntax changed MAY-2014
    error::warnAboutAge("constraint keyword", 1406);

    return usable;
}


void Foam::decompositionMethod::readConstraints()
{
    constraints_.clear();

    const dictionary* dictptr = decompDict_.findDict("constraints");

    if (dictptr)
    {
        for (const entry& dEntry : *dictptr)
        {
            if (!dEntry.isDict())  // safety
            {
                // Ignore or warn
                continue;
            }

            const dictionary& dict = dEntry.dict();

            if (dict.getOrDefault("enabled", true))
            {
                constraints_.append(decompositionConstraint::New(dict));
            }
        }
    }

    // Backwards compatibility (MAY-2014)
    if (constraintCompat("preserveBaffles"))
    {
        constraints_.append
        (
            new decompositionConstraints::preserveBaffles()
        );
    }

    if (constraintCompat("preservePatches"))
    {
        constraints_.append
        (
            new decompositionConstraints::preservePatches
            (
                decompDict_.get<wordRes>("preservePatches")
            )
        );
    }

    if (constraintCompat("preserveFaceZones"))
    {
        constraints_.append
        (
            new decompositionConstraints::preserveFaceZones
            (
                decompDict_.get<wordRes>("preserveFaceZones")
            )
        );
    }

    if (constraintCompat("singleProcessorFaceSets"))
    {
        constraints_.append
        (
            new decompositionConstraints::singleProcessorFaceSets
            (
                decompDict_.lookup("singleProcessorFaceSets")
            )
        );
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::decompositionMethod::findCoeffsDict
(
    const dictionary& dict,
    const word& coeffsName,
    int select
)
{
    const bool allowDefault = !(select & selectionType::EXACT);

    const dictionary* dictptr =
        cfindCoeffsDict(dict, coeffsName, allowDefault);

    if (dictptr)
    {
        return *dictptr;
    }

    // Not found
    if (select & selectionType::MANDATORY)
    {
        FatalIOError
            << "'" << coeffsName << "' dictionary not found in dictionary "
            << dict.name() << endl
            << abort(FatalIOError);
    }

    if (select & selectionType::NULL_DICT)
    {
        return dictionary::null;
    }

    return dict;  // Return original dictionary
}


const Foam::dictionary& Foam::decompositionMethod::findCoeffsDict
(
    const word& coeffsName,
    int select
) const
{
    const bool allowDefault = !(select & selectionType::EXACT);

    const dictionary* dictptr = nullptr;

    if (!decompRegionDict_.empty())
    {
        // Region-specific dictionary
        dictptr = cfindCoeffsDict(decompRegionDict_, coeffsName, allowDefault);
    }
    if (!dictptr)
    {
        // General
        dictptr = cfindCoeffsDict(decompDict_, coeffsName, allowDefault);
    }

    if (dictptr)
    {
        return *dictptr;
    }

    // Not found
    if (select & selectionType::MANDATORY)
    {
        FatalIOError
            << "'" << coeffsName << "' dictionary not found in dictionary "
            << decompDict_.name() << endl
            << abort(FatalIOError);
    }

    if (select & selectionType::NULL_DICT)
    {
        return dictionary::null;
    }

    return decompDict_;  // Return general dictionary
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethod::decompositionMethod(const label numDomains)
:
    decompDict_(dictionary::null),
    decompRegionDict_(dictionary::null),
    nDomains_(numDomains)
{}


Foam::decompositionMethod::decompositionMethod
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompDict_(decompDict),
    decompRegionDict_
    (
        optionalRegionDict(decompDict_, regionName)
    ),
    nDomains_(nDomains(decompDict, regionName))
{
    readConstraints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompDict,
    const word& regionName
)
{
    word methodType(decompDict.get<word>("method"));

    const dictionary& regionDict = optionalRegionDict(decompDict, regionName);
    regionDict.readIfPresent("method", methodType);

    auto* ctorPtr = dictionaryConstructorTable(methodType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            decompDict,
            "decompositionMethod",
            methodType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // verbose
    {
        Info<< "Decomposition method " << methodType
            << " [" << (nDomains(decompDict, regionName)) << ']';

        if (!regionName.empty())
        {
            Info<< " (region " << regionName << ')';
        }
        Info<< endl;
    }

    return autoPtr<decompositionMethod>(ctorPtr(decompDict, regionName));
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
) const
{
    CompactListList<label> coarseCellCells;
    globalMeshData::calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),
        true,        // Global mesh connectivity
        coarseCellCells
    );

    // Decompose based on agglomerated points
    labelList decomp
    (
        decompose
        (
            coarseCellCells,
            coarsePoints,
            coarseWeights
        )
    );

    // From coarse back to fine for original mesh
    return labelList(decomp, fineToCoarse);
}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells
)
{
    globalMeshData::calcCellCells
    (
        mesh,
        agglom,
        nLocalCoarse,
        parallel,
        cellCells
    );
}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells,
    CompactListList<scalar>& cellCellWeights
)
{
    globalMeshData::calcCellCells
    (
        mesh,
        agglom,
        nLocalCoarse,
        parallel,
        cellCells,
        cellCellWeights
    );
}


// NOTE:
// - alternative calcCellCells that handled explicitConnections was
//   deactivated (2014 or earlier) and finally removed APR-2018.

Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const scalarField& cellWeights,

    //- Whether owner and neighbour should be on same processor
    //  (takes priority over explicitConnections)
    const boolList& blockedFace,

    //- Whether whole sets of faces (and point neighbours) need to be kept
    //  on single processor
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,

    //- Additional connections between boundary faces
    const List<labelPair>& explicitConnections
) const
{
    // Any weights specified?
    const bool hasWeights = returnReduceOr(cellWeights.size());

    if (hasWeights && (cellWeights.size() != mesh.nCells()))
    {
        FatalErrorInFunction
            << "Number of weights (" << cellWeights.size()
            << ") != number of cells (" << mesh.nCells() << ")"
            << exit(FatalError);
    }

    // Any faces not blocked?
    const bool hasUnblocked =
        returnReduceOr
        (
            !blockedFace.empty() && !BitOps::all(blockedFace)
        );


    // Any non-mesh connections?
    const label nConnections = returnReduce
    (
        explicitConnections.size(),
        sumOp<label>()
    );


    // Any processor sets?
    label nProcSets = 0;
    for (const labelList& procset : specifiedProcessorFaces)
    {
        nProcSets += procset.size();
    }
    reduce(nProcSets, sumOp<label>());


    // Either do decomposition on cell centres or on agglomeration

    if (!hasUnblocked && !nConnections && !nProcSets)
    {
        // No constraints, possibly weights

        return
        (
            hasWeights
          ? decompose(mesh, mesh.cellCentres(), cellWeights)
          : decompose(mesh, mesh.cellCentres())
        );
    }


    // The harder work.
    // When we have processor sets, connections, or blocked faces.


    // Determine local regions, separated by blockedFaces
    regionSplit localRegion(mesh, blockedFace, explicitConnections, false);

    if (debug)
    {
        // Only need to count unblocked faces for debugging
        const label nUnblocked =
        (
            hasUnblocked
          ? returnReduce
            (
                label(BitOps::count(blockedFace, false)),
                sumOp<label>()
            )
          : 0
        );

        Info<< "Constrained decomposition:" << nl
            << "    faces with same owner and neighbour processor : "
            << nUnblocked << nl
            << "    baffle faces with same owner processor        : "
            << nConnections << nl
            << "    faces all on same processor                   : "
            << nProcSets << nl
            << "    split into " << localRegion.nLocalRegions()
            << " regions."
            << endl;
    }


    // Gather region weights and determine region cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // For the region centre, just take the first cell in the region.
    // If we average the region centre instead, cyclics could cause
    // the average domain centre to be outside of domain.

    scalarField regionWeights(localRegion.nLocalRegions(), Foam::zero{});

    pointField regionCentres(localRegion.nLocalRegions(), point::max);

    if (hasWeights)
    {
        forAll(localRegion, celli)
        {
            const label regioni = localRegion[celli];

            regionWeights[regioni] += cellWeights[celli];

            if (regionCentres[regioni] == point::max)
            {
                regionCentres[regioni] = mesh.cellCentres()[celli];
            }
        }
    }
    else
    {
        forAll(localRegion, celli)
        {
            const label regioni = localRegion[celli];

            regionWeights[regioni] += 1.0;

            if (regionCentres[regioni] == point::max)
            {
                regionCentres[regioni] = mesh.cellCentres()[celli];
            }
        }
    }

    // Do decomposition on agglomeration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList finalDecomp =
        decompose
        (
            mesh,
            localRegion,
            regionCentres,
            regionWeights
        );


    // Apply explicitConnections since decompose did not know about them
    for (const labelPair& baffle : explicitConnections)
    {
        const label f0 = baffle.first();
        const label f1 = baffle.second();

        if (!blockedFace[f0] && !blockedFace[f1])
        {
            // Note: what if internal faces and owner and neighbour on
            // different processor?
            // So for now just push owner side proc

            const label proci = finalDecomp[mesh.faceOwner()[f0]];

            finalDecomp[mesh.faceOwner()[f1]] = proci;
            if (mesh.isInternalFace(f1))
            {
                finalDecomp[mesh.faceNeighbour()[f1]] = proci;
            }
        }
        else if (blockedFace[f0] != blockedFace[f1])
        {
            FatalErrorInFunction
                << "On explicit connection between faces " << f0
                << " and " << f1
                << " the two blockedFace status are not equal : "
                << blockedFace[f0] << " and " << blockedFace[f1]
                << exit(FatalError);
        }
    }


    // blockedFaces corresponding to processor faces need to be handled
    // separately since not handled by local regionSplit. We need to
    // walk now across coupled faces and make sure to move a whole
    // global region across

    // This additionally consolidates/compacts the regions numbers globally,
    // since that was skipped in the previous regionSplit.
    if (Pstream::parRun())
    {
        // Re-do regionSplit

        // Field on cells and faces.
        List<minData> cellData(mesh.nCells());
        List<minData> faceData(mesh.nFaces());

        // Take over blockedFaces by seeding a negative number
        // (so is always less than the decomposition)
        label nUnblocked = 0;
        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                faceData[facei] = minData(-123);
            }
            else
            {
                ++nUnblocked;
            }
        }

        // Seed unblocked faces with destination processor
        labelList seedFaces(nUnblocked);
        List<minData> seedData(nUnblocked);
        nUnblocked = 0;

        forAll(blockedFace, facei)
        {
            if (!blockedFace[facei])
            {
                const label own = mesh.faceOwner()[facei];
                seedFaces[nUnblocked] = facei;
                seedData[nUnblocked] = minData(finalDecomp[own]);
                nUnblocked++;
            }
        }


        // Propagate information inwards
        FaceCellWave<minData> deltaCalc
        (
            mesh,
            seedFaces,
            seedData,
            faceData,
            cellData,
            mesh.globalData().nTotalCells()+1
        );

        // And extract
        forAll(finalDecomp, celli)
        {
            if (cellData[celli].valid(deltaCalc.data()))
            {
                finalDecomp[celli] = cellData[celli].data();
            }
        }
    }


    // For specifiedProcessorFaces rework the cellToProc to enforce
    // all on one processor since we can't guarantee that the input
    // to regionSplit was a single region.
    // E.g. faceSet 'a' with the cells split into two regions
    // by a notch formed by two walls
    //
    //          \   /
    //           \ /
    //    ---a----+-----a-----
    //
    //
    // Note that reworking the cellToProc might make the decomposition
    // unbalanced.
    forAll(specifiedProcessorFaces, seti)
    {
        const labelList& set = specifiedProcessorFaces[seti];

        label proci = specifiedProcessor[seti];
        if (proci == -1)
        {
            // If no processor specified - use the one from the 0th element
            if (set.size())
            {
                proci = finalDecomp[mesh.faceOwner()[set[0]]];
            }
            else
            {
                // Zero-sized processor (e.g. from redistributePar)
                proci = 0;
            }
        }

        for (const label facei : set)
        {
            const face& f = mesh.faces()[facei];
            for (const label pointi : f)
            {
                const labelList& pFaces = mesh.pointFaces()[pointi];
                for (const label pFacei : pFaces)
                {
                    finalDecomp[mesh.faceOwner()[pFacei]] = proci;
                    if (mesh.isInternalFace(pFacei))
                    {
                        finalDecomp[mesh.faceNeighbour()[pFacei]] = proci;
                    }
                }
            }
        }
    }


    if (debug && Pstream::parRun())
    {
        labelList nbrDecomp;
        syncTools::swapBoundaryCellList(mesh, finalDecomp, nbrDecomp);

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    const label facei = pp.start()+i;
                    const label own = mesh.faceOwner()[facei];
                    const label bFacei = facei-mesh.nInternalFaces();

                    if (!blockedFace[facei])
                    {
                        const label ownProc = finalDecomp[own];
                        const label nbrProc = nbrDecomp[bFacei];

                        if (ownProc != nbrProc)
                        {
                            FatalErrorInFunction
                                << "patch:" << pp.name()
                                << " face:" << facei
                                << " at:" << mesh.faceCentres()[facei]
                                << " ownProc:" << ownProc
                                << " nbrProc:" << nbrProc
                                << exit(FatalError);
                        }
                    }
                }
            }
        }
    }

    return finalDecomp;
}


void Foam::decompositionMethod::setConstraints
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.resize_nocopy(mesh.nFaces());
    blockedFace = true;

    specifiedProcessorFaces.clear();
    explicitConnections.clear();

    for (const decompositionConstraint& decompConstraint : constraints_)
    {
        decompConstraint.add
        (
            mesh,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections
        );
    }
}


void Foam::decompositionMethod::applyConstraints
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    for (const decompositionConstraint& decompConstraint : constraints_)
    {
        decompConstraint.apply
        (
            mesh,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections,
            decomposition
        );
    }
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const scalarField& cellWeights
) const
{
    // Collect all constraints

    boolList blockedFace;
    PtrList<labelList> specifiedProcessorFaces;
    labelList specifiedProcessor;
    List<labelPair> explicitConnections;
    setConstraints
    (
        mesh,
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections
    );


    // Construct decomposition method and either do decomposition on
    // cell centres or on agglomeration

    labelList finalDecomp = decompose
    (
        mesh,
        cellWeights,            // optional weights
        blockedFace,            // any cells to be combined
        specifiedProcessorFaces,// any whole cluster of cells to be kept
        specifiedProcessor,
        explicitConnections     // baffles
    );


    // Give any constraint the option of modifying the decomposition

    applyConstraints
    (
        mesh,
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections,
        finalDecomp
    );

    return finalDecomp;
}


// * * * * * * * * * * * * * * * Stub Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethod::decompose
(
    const pointField& points,
    const scalarField& pointWeights
) const
{
    NotImplemented;
    return labelList();
}


// ************************************************************************* //
