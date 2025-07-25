/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Wikki Ltd
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "faMesh.H"
#include "faPatchData.H"
#include "emptyFaPatch.H"
#include "ignoreFaPatch.H"
#include "processorFaPatch.H"
#include "processorPolyPatch.H"
#include "foamVtkLineWriter.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Write edges in VTK format
template<class PatchType>
static void vtkWritePatchEdges
(
    const PatchType& p,
    const labelList& selectEdges,
    const fileName& outputPath,
    const word& outputName
)
{
    edgeList dumpEdges(p.edges(), selectEdges);

    vtk::lineWriter writer
    (
        p.localPoints(),
        dumpEdges,
        outputPath/outputName
    );

    writer.writeGeometry();

    // CellData
    writer.beginCellData();
    writer.writeProcIDs();
    writer.close();
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMesh::addFaPatches
(
    faPatchList& plist,
    const bool validBoundary
)
{
    if (!boundary().empty())
    {
        FatalErrorInFunction
            << "boundary already exists"
            << abort(FatalError);
    }

    globalMeshDataPtr_.reset(nullptr);

    boundary_.transfer(plist);

    setPrimitiveMeshData();

    if (validBoundary)
    {
        boundary_.checkDefinition();
    }
}


void Foam::faMesh::addFaPatches
(
    const List<faPatch*>& p,
    const bool validBoundary
)
{
    // Acquire ownership of the pointers
    faPatchList plist(const_cast<List<faPatch*>&>(p));

    addFaPatches(plist, validBoundary);
}


Foam::faPatchList Foam::faMesh::createOnePatch
(
    const word& patchName,
    const word& patchType
) const
{
    dictionary onePatchDict;
    if (!patchName.empty())
    {
        onePatchDict.add("name", patchName);
    }
    if (!patchType.empty())
    {
        onePatchDict.add("type", patchType);
    }

    return createPatchList
    (
        dictionary::null,
        word::null,     // Name for empty patch placeholder
        &onePatchDict   // Definitions for defaultPatch
    );
}


Foam::faPatchList Foam::faMesh::createPatchList
(
    const dictionary& bndDict,
    const word& emptyPatchName,
    const dictionary* defaultPatchDefinition
) const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    // Transcribe into patch definitions
    DynamicList<faPatchData> faPatchDefs(bndDict.size() + 8);
    for (const entry& dEntry : bndDict)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Not a dictionary entry: " << dEntry.name() << nl;
            continue;
        }
        const dictionary& patchDict = dEntry.dict();

        // Add entry
        auto& patchDef = faPatchDefs.emplace_back();
        patchDef.name_ = dEntry.keyword();
        patchDef.type_ = patchDict.get<word>("type");

        word patchName;

        // Optional: ownerPolyPatch
        if (patchDict.readIfPresent("ownerPolyPatch", patchName))
        {
            patchDef.ownerPolyPatchId_ = pbm.findPatchID(patchName);
            if (patchDef.ownerPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "ownerPolyPatch " << patchName << " not found"
                    << exit(FatalError);
            }
        }

        // Mandatory: neighbourPolyPatch
        patchDict.readEntry("neighbourPolyPatch", patchName);
        {
            patchDef.neighPolyPatchId_ = pbm.findPatchID(patchName);
            if (patchDef.neighPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "neighbourPolyPatch " << patchName << " not found"
                    << exit(FatalError);
            }
        }
    }

    // Additional empty placeholder patch?
    if (!emptyPatchName.empty())
    {
        auto& patchDef = faPatchDefs.emplace_back();
        patchDef.name_ = emptyPatchName;
        patchDef.type_ = emptyFaPatch::typeName_();
    }

    label nWarnUndefinedPatch(5);

    // Placeholder for any undefined edges
    const label undefPatchIndex = faPatchDefs.size();
    {
        auto& patchDef = faPatchDefs.emplace_back();
        patchDef.name_ = "undefined";
        patchDef.type_ = "patch";

        if (defaultPatchDefinition)
        {
            if
            (
                (*defaultPatchDefinition).readIfPresent("name", patchDef.name_)
            )
            {
                // Suppress warnings if defaultPatch name was specified
                // - probably means we want to use this mechanism
                nWarnUndefinedPatch = 0;
            }
            (*defaultPatchDefinition).readIfPresent("type", patchDef.type_);
        }
    }

    // Placeholder for any undefined edges
    const label ignorePatchIndex = faPatchDefs.size();
    {
        auto& patchDef = faPatchDefs.emplace_back();
        patchDef.name_ = "_ignore_edges_";
        patchDef.type_ = ignoreFaPatch::typeName_();
    }

    // ----------------------------------------------------------------------

    // Get edge connections (in globally consistent ordering)
    List<Pair<patchTuple>> bndEdgeConnections
    (
        this->getBoundaryEdgeConnections()
    );

    // Update accordingly
    this->setBoundaryConnections(bndEdgeConnections);


    // Lookup of patchDef for each connection. Initially all -1 (unvisited)
    labelList patchDefLookup(bndEdgeConnections.size(), -1);

    Map<labelHashSet> procConnections;
    labelHashSet patchDefsUsed;

    labelHashSet badEdges(2*bndEdgeConnections.size());

    forAll(bndEdgeConnections, connecti)
    {
        const Pair<patchTuple>& connection = bndEdgeConnections[connecti];
        const auto& a = connection.first();
        const auto& b = connection.second();

        edge patchPair;

        if (!a.valid() || !b.valid())
        {
            // If either is invalid, mark as an 'ignore' edge
            patchDefLookup[connecti] = ignorePatchIndex;
            patchDefsUsed.insert(ignorePatchIndex);
            continue;
        }
        else if (a.is_finiteArea())
        {
            if (b.is_finiteArea())
            {
                // Expecting an inter-processor connection

                if (a.procNo() == b.procNo())
                {
                    // An intra-processor connection (should not be possible)
                    FatalErrorInFunction
                        << "Processor-processor addressing error:" << nl
                        << "Both connections have the same processor: "
                        << a.procNo() << nl
                        << "Connecting patches "
                        << a.realPatchi() << " and " << b.realPatchi() << nl
                        << abort(FatalError);
                }
                else if (a.is_localProc())
                {
                    procConnections(b.procNo()).insert(connecti);
                }
                else
                {
                    procConnections(a.procNo()).insert(connecti);
                }

                continue;
            }
            else if (a.is_localProc())
            {
                patchPair.first()  = a.realPatchi();
                patchPair.second() = b.realPatchi();
            }
        }
        else if (b.is_finiteArea() && b.is_localProc())
        {
            patchPair.first()  = b.realPatchi();
            patchPair.second() = a.realPatchi();
        }


        // Find first definition with a matching neighbour and
        // possibly with a matching owner.

        label bestPatchDefi = -1;

        const label nPatchDefs = (patchPair.valid() ? faPatchDefs.size() : 0);

        for (label patchDefi = 0; patchDefi < nPatchDefs; ++patchDefi)
        {
            const int match = faPatchDefs[patchDefi].matchPatchPair(patchPair);
            if (match == 3)
            {
                // Exact match (owner/neighbour) - done!
                bestPatchDefi = patchDefi;
                break;
            }
            else if (match == 2 && bestPatchDefi < 0)
            {
                // Match (neighbour) - keep looking for exact match
                bestPatchDefi = patchDefi;
            }
        }

        if (bestPatchDefi < 0)
        {
            bestPatchDefi = undefPatchIndex;  // Missed?
        }

        patchDefLookup[connecti] = bestPatchDefi;
        patchDefsUsed.insert(bestPatchDefi);
    }


    bool reportBadEdges = false;

    // Skip undefPatchIndex if not actually needed anywhere
    if (!returnReduceOr(patchDefsUsed.found(undefPatchIndex)))
    {
        faPatchDefs[undefPatchIndex].clear();
    }
    else
    {
        patchDefsUsed.insert(undefPatchIndex);  // Parallel consistency
        reportBadEdges = true;
    }

    // Skip ignorePatchIndex if not actually needed anywhere
    if (!returnReduceOr(patchDefsUsed.found(ignorePatchIndex)))
    {
        faPatchDefs[ignorePatchIndex].clear();
    }
    else
    {
        patchDefsUsed.insert(ignorePatchIndex);  // Parallel consistency
        reportBadEdges = true;
    }

    // Report locations of undefined edges
    if (reportBadEdges)
    {
        badEdges.clear();
        forAll(patchDefLookup, connecti)
        {
            if (patchDefLookup[connecti] == undefPatchIndex)
            {
                const auto& connection = bndEdgeConnections[connecti];

                const auto& a = connection.first();
                const auto& b = connection.second();

                if (a.is_localProc() && a.is_finiteArea())
                {
                    badEdges.insert(a.patchEdgei());
                }
                else if (b.is_localProc() && b.is_finiteArea())
                {
                    badEdges.insert(b.patchEdgei());
                }

                if (badEdges.size() <= nWarnUndefinedPatch)
                {
                    Pout<< "Undefined connection: "
                        << "(patch:" << a.realPatchi()
                        << " face:" << a.meshFacei()
                        << " proc:" << a.procNo()
                        << ") and (patch:" << b.realPatchi()
                        << " face:" << b.meshFacei()
                        << " proc:" << b.procNo()
                        << ")  patch:"
                        <<
                        (
                            a.realPatchi() >= 0
                          ? pbm[a.realPatchi()].name()
                          : word::null
                        )
                        << " and patch:"
                        <<
                        (
                            b.realPatchi() >= 0
                          ? pbm[b.realPatchi()].name()
                          : word::null
                        )
                        << nl;
                }
            }
        }

        if (returnReduceOr(badEdges.size()))
        {
            // Report directly as Info, not InfoInFunction
            // since it can also be an expected result when
            // nWarnUndefinedPatch == 0
            Info<< nl
                << "Had "
                << returnReduce(badEdges.size(), sumOp<label>()) << '/'
                << returnReduce(patch().nBoundaryEdges(), sumOp<label>())
                << " undefined edge connections, added to defaultPatch: "
                << faPatchDefs[undefPatchIndex].name_ << nl << nl
                << "==> Could indicate a non-manifold patch geometry" << nl
                << nl;

            if (nWarnUndefinedPatch)
            {
                labelList selectEdges(badEdges.sortedToc());
                word outputName("faMesh-construct.undefEdges");

                vtkWritePatchEdges
                (
                    patch(),
                    selectEdges,
                    thisDb().time().globalPath(),
                    outputName
                );

                InfoInFunction
                    << "(debug) wrote " << outputName << nl;
            }
        }
    }

    // Report locations of undefined edges
    if (reportBadEdges)
    {
        badEdges.clear();
        forAll(patchDefLookup, connecti)
        {
            if (patchDefLookup[connecti] == ignorePatchIndex)
            {
                const auto& connection = bndEdgeConnections[connecti];

                const auto& a = connection.first();
                const auto& b = connection.second();

                if (a.is_localProc() && a.is_finiteArea())
                {
                    badEdges.insert(a.patchEdgei());
                }
                else if (b.is_localProc() && b.is_finiteArea())
                {
                    badEdges.insert(b.patchEdgei());
                }

                if (badEdges.size() <= nWarnUndefinedPatch)
                {
                    Pout<< "Illegal connection: "
                        << "(patch:" << a.realPatchi()
                        << " face:" << a.meshFacei()
                        << " proc:" << a.procNo()
                        << ") and (patch:" << b.realPatchi()
                        << " face:" << b.meshFacei()
                        << " proc:" << b.procNo()
                        << ") patch:"
                        <<
                        (
                            a.realPatchi() >= 0
                          ? pbm[a.realPatchi()].name()
                          : word::null
                        )
                        << " and patch:"
                        <<
                        (
                            b.realPatchi() >= 0
                          ? pbm[b.realPatchi()].name()
                          : word::null
                        )
                        << nl;
                }
            }
        }

        if (returnReduceOr(badEdges.size()))
        {
            // Report directly as Info, not InfoInFunction
            // since it can also be an expected result when
            // nWarnUndefinedPatch == 0
            Info<< nl
                << "Had "
                << returnReduce(badEdges.size(), sumOp<label>()) << '/'
                << returnReduce(patch().nBoundaryEdges(), sumOp<label>())
                << " illegal edge connections, added to "
                << faPatchDefs[ignorePatchIndex].name_ << nl << nl
                << "==> Could indicate a non-manifold patch geometry" << nl
                << nl;

            if (nWarnUndefinedPatch)
            {
                labelList selectEdges(badEdges.sortedToc());
                word outputName("faMesh-construct.ignoreEdges");

                vtkWritePatchEdges
                (
                    patch(),
                    selectEdges,
                    thisDb().time().globalPath(),
                    outputName
                );

                InfoInFunction
                    << "(debug) wrote " << outputName << nl;
            }
        }
    }

    // Create processor-processor definitions
    Map<label> procToDefLookup;
    {
        procToDefLookup.reserve(procConnections.size());
        faPatchDefs.reserve(faPatchDefs.size() + procConnections.size());

        for (const label otherProci : procConnections.sortedToc())
        {
            const label patchDefi = faPatchDefs.size();
            procToDefLookup.insert(otherProci, patchDefi);

            // Add entry
            auto& patchDef = faPatchDefs.emplace_back();

            patchDef.assign_coupled(UPstream::myProcNo(), otherProci);
        }
    }


    // Extract which edges map to which patch
    // and set edgeLabels for each faPatch

    DynamicList<label> selectEdges(bndEdgeConnections.size());
    label nOffProcessorEdges = 0;

    for (const label patchDefi : patchDefsUsed.sortedToc())
    {
        auto& patchDef = faPatchDefs[patchDefi];
        selectEdges.clear();

        // Find the corresponding entries
        // and extract the patchEdgeId

        forAll(patchDefLookup, connecti)
        {
            if (patchDefLookup[connecti] == patchDefi)
            {
                const auto& a = bndEdgeConnections[connecti].first();
                const auto& b = bndEdgeConnections[connecti].second();

                if (a.is_localProc() && a.is_finiteArea())
                {
                    selectEdges.push_back(a.patchEdgei());
                }
                else if (b.is_localProc() && b.is_finiteArea())
                {
                    selectEdges.push_back(b.patchEdgei());
                }
                else if (a.valid() && b.valid())
                {
                    FatalErrorInFunction
                        << "Error in programming logic" << nl
                        << abort(FatalError);
                }

                if (a.is_localProc() != b.is_localProc())
                {
                    ++nOffProcessorEdges;
                }

                patchDefLookup[connecti] = -2;  // Mark as handled
            }
        }

        // Remove any cosmetic sorting artifacts from off-processor
        // termination by doing using a regular sort here.

        Foam::sort(selectEdges);
        patchDef.edgeLabels_ = selectEdges;
    }

    if (debug)
    {
        Pout<< "Had " << nOffProcessorEdges
            << " patch edges connected off-processor" << endl;

        InfoInFunction
            << "Total "
            << returnReduce(nOffProcessorEdges, sumOp<label>())
            << " patch edges connected off-processor" << endl;
    }


    // Processor patches
    for (const label otherProci : procToDefLookup.sortedToc())
    {
        const label patchDefi = procToDefLookup[otherProci];

        auto& patchDef = faPatchDefs[patchDefi];
        selectEdges.clear();

        // Find the corresponding entries
        // and extract the patchEdgeId

        for (const label connecti : procConnections(otherProci).sortedToc())
        {
            const auto& connection = bndEdgeConnections[connecti];
            const auto& a = connection.first();
            const auto& b = connection.second();

            if (a.is_localProc())
            {
                selectEdges.push_back(a.patchEdgei());
            }
            else if (b.is_localProc())
            {
                selectEdges.push_back(b.patchEdgei());
            }
            else
            {
                FatalErrorInFunction
                    << "Error in programming logic" << nl
                    << abort(FatalError);
            }

            patchDefLookup[connecti] = -2;  // Mark as handled
        }

        // The edge order is guaranteed to be consistent from the original
        // getBoundaryEdgeConnections() - sorted by proc/patch/edge

        patchDef.edgeLabels_ = selectEdges;
    }


    // Now convert list of definitions to list of patches

    label nPatches = 0;
    faPatchList newPatches(faPatchDefs.size());

    for (faPatchData& patchDef : faPatchDefs)
    {
        if (!patchDef.good())
        {
            continue;
        }

        newPatches.set
        (
            nPatches,
            faPatch::New
            (
                patchDef.name(),        // name
                patchDef.dict(false),   // withEdgeLabels == false
                nPatches,               // index
                boundary()
            )
        );

        // Transfer edge labels
        newPatches[nPatches].resetEdges(std::move(patchDef.edgeLabels_));
        ++nPatches;
    }
    newPatches.resize(nPatches);

    if (debug > 1)
    {
        label nTotal = 0;
        Pout<< "Created new finiteArea patches:" << nl;
        for (const faPatch& p : newPatches)
        {
            Pout<< "  size: " << p.size()
                << " name:" << p.name()
                << " type:" << p.type() << nl;
            nTotal += p.size();
        }

        Pout<< "addressed: " << nTotal
            << '/' << patch().nBoundaryEdges() << " edges" << endl;
    }

    return newPatches;
}


// ************************************************************************* //
