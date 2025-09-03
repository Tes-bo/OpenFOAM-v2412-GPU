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

#include "refinementHistory.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinementHistory, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::refinementHistory::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


void Foam::refinementHistory::writeEntry
(
    const List<splitCell8>& splitCells,
    const splitCell8& split
)
{
    // Write me:
    if (split.addedCellsPtr_)
    {
        Pout<< "parent:" << split.parent_
            << " subCells:" << split.addedCellsPtr_()
            << endl;
    }
    else
    {
        Pout<< "parent:" << split.parent_
            << " no subcells"
            << endl;
    }

    if (split.parent_ >= 0)
    {
        Pout<< "parent data:" << endl;
        // Write my parent
        string oldPrefix = Pout.prefix();
        Pout.prefix() = "  " + oldPrefix;
        writeEntry(splitCells, splitCells[split.parent_]);
        Pout.prefix() = oldPrefix;
    }
}


void Foam::refinementHistory::writeDebug
(
    const labelList& visibleCells,
    const List<splitCell8>& splitCells
)
{
    string oldPrefix = Pout.prefix();
    Pout.prefix() = "";

    forAll(visibleCells, celli)
    {
        label index = visibleCells[celli];

        if (index >= 0)
        {
            Pout<< "Cell from refinement:" << celli << " index:" << index
                << endl;

            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[index]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            Pout<< "Unrefined cell:" << celli << " index:" << index << endl;
        }
    }
    Pout.prefix() = oldPrefix;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistory::splitCell8::splitCell8()
:
    parent_(-1),
    addedCellsPtr_(nullptr)
{}


Foam::refinementHistory::splitCell8::splitCell8(const label parent)
:
    parent_(parent),
    addedCellsPtr_(nullptr)
{}


Foam::refinementHistory::splitCell8::splitCell8(Istream& is)
:
    splitCell8()
{
    is >> *this;
}


Foam::refinementHistory::splitCell8::splitCell8(const splitCell8& rhs)
:
    parent_(rhs.parent_),
    addedCellsPtr_(rhs.addedCellsPtr_.clone()) // Copy, not steal
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::refinementHistory::splitCell8::operator=(const splitCell8& rhs)
{
    // Assignment operator since autoPtr otherwise 'steals' storage.

    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    parent_ = rhs.parent_;
    addedCellsPtr_.reset(rhs.addedCellsPtr_.clone()); // Copy, not steal
}


bool Foam::refinementHistory::splitCell8::operator==
(
    const splitCell8& rhs
)
const
{
    if (parent_ != rhs.parent_)
    {
        return false;
    }
    if (bool(addedCellsPtr_) != bool(rhs.addedCellsPtr_))
    {
        return false;
    }
    else if (addedCellsPtr_)  // With previous test, means rhs is also defined
    {
        return addedCellsPtr_() == rhs.addedCellsPtr_();
    }

    return true;
}


bool Foam::refinementHistory::splitCell8::operator!=
(
    const splitCell8& rhs
) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistory::splitCell8& sc)
{
    labelList addedCells;

    is >> sc.parent_ >> addedCells;

    if (addedCells.size())
    {
        sc.addedCellsPtr_.reset(new FixedList<label, 8>(addedCells));
    }
    else
    {
        sc.addedCellsPtr_.reset(nullptr);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementHistory::splitCell8& sc
)
{
    // Output as labelList so we can have 0 sized lists. Alternative is to
    // output as fixedlist with e.g. -1 elements and check for this upon
    // reading. However would cause much more data to be transferred.

    labelList labels;

    if (sc.addedCellsPtr_)
    {
        labels = sc.addedCellsPtr_();
    }

    return os << sc.parent_ << token::SPACE << labels;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementHistory::checkIndices() const
{
    // Check indices.
    forAll(visibleCells_, i)
    {
        if (visibleCells_[i] < 0 && visibleCells_[i] >= splitCells_.size())
        {
            FatalErrorInFunction
                << "Illegal entry " << visibleCells_[i]
                << " in visibleCells at location" << i << nl
                << "It points outside the range of splitCells : 0.."
                << splitCells_.size()-1
                << abort(FatalError);
        }
    }
}


Foam::label Foam::refinementHistory::allocateSplitCell
(
    const label parent,
    const label i
)
{
    label index = -1;

    if (freeSplitCells_.size())
    {
        index = freeSplitCells_.back();
        freeSplitCells_.pop_back();

        splitCells_[index] = splitCell8(parent);
    }
    else
    {
        index = splitCells_.size();
        splitCells_.push_back(splitCell8(parent));
    }


    // Update the parent field
    if (parent >= 0)
    {
        splitCell8& parentSplit = splitCells_[parent];

        if (!parentSplit.addedCellsPtr_)
        {
            // Allocate storage on parent for the 8 subcells.
            parentSplit.addedCellsPtr_.reset(new FixedList<label, 8>(-1));
        }


        // Store me on my parent
        FixedList<label, 8>& parentSplits = parentSplit.addedCellsPtr_();

        parentSplits[i] = index;
    }

    return index;
}


void Foam::refinementHistory::freeSplitCell(const label index)
{
    splitCell8& split = splitCells_[index];

    // Make sure parent does not point to me anymore.
    if (split.parent_ >= 0)
    {
        autoPtr<FixedList<label, 8>>& subCellsPtr =
            splitCells_[split.parent_].addedCellsPtr_;

        if (subCellsPtr)
        {
            FixedList<label, 8>& subCells = subCellsPtr();

            label myPos = subCells.find(index);

            if (myPos == -1)
            {
                FatalErrorInFunction
                    << "Problem: cannot find myself in"
                    << " parents' children" << abort(FatalError);
            }
            else
            {
                subCells[myPos] = -1;
            }
        }
    }

    // Mark splitCell as free
    split.parent_ = -2;

    // Add to cache of free splitCells
    freeSplitCells_.append(index);
}


void Foam::refinementHistory::markSplit
(
    const label index,
    labelList& oldToNew,
    DynamicList<splitCell8>& newSplitCells
) const
{
    if (oldToNew[index] == -1)
    {
        // Not yet compacted.

        const splitCell8& split = splitCells_[index];

        oldToNew[index] = newSplitCells.size();
        newSplitCells.append(split);

        if (split.parent_ >= 0)
        {
            markSplit(split.parent_, oldToNew, newSplitCells);
        }
        if (split.addedCellsPtr_)
        {
            const FixedList<label, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    markSplit(splits[i], oldToNew, newSplitCells);
                }
            }
        }
    }
}


void Foam::refinementHistory::mark
(
    const label val,
    const label index,
    labelList& splitToVal
) const
{
    splitToVal[index] = val;

    const splitCell8& split = splitCells_[index];

    if (split.addedCellsPtr_)
    {
        const FixedList<label, 8>& splits = split.addedCellsPtr_();

        forAll(splits, i)
        {
            if (splits[i] >= 0)
            {
                mark(val, splits[i], splitToVal);
            }
        }
    }
}


Foam::label Foam::refinementHistory::markCommonCells
(
    labelList& cellToCluster
) const
{
    label clusterI = 0;

    labelList splitToCluster(splitCells_.size(), -1);

    // Pass1: find top of all clusters
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Find highest ancestor
            while (splitCells_[index].parent_ != -1)
            {
                index = splitCells_[index].parent_;
            }

            // Mark tree with clusterI
            if (splitToCluster[index] == -1)
            {
                mark(clusterI, index, splitToCluster);
                clusterI++;
            }
        }
    }

    // Pass2: mark all cells with cluster
    cellToCluster.setSize(visibleCells_.size(), -1);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            cellToCluster[cellI] = splitToCluster[index];
        }
    }

    return clusterI;
}


void Foam::refinementHistory::add
(
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    blockedFace.setSize(mesh.nFaces(), true);

    // Find common parent for all cells
    labelList cellToCluster;
    markCommonCells(cellToCluster);


    // Unblock all faces inbetween same cluster

    label nUnblocked = 0;

    forAll(mesh.faceNeighbour(), faceI)
    {
        label ownCluster = cellToCluster[mesh.faceOwner()[faceI]];
        label neiCluster = cellToCluster[mesh.faceNeighbour()[faceI]];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (blockedFace[faceI])
            {
                blockedFace[faceI] = false;
                nUnblocked++;
            }
        }
    }

    if (refinementHistory::debug)
    {
        Info<< type() << " : unblocked "
            << returnReduce(nUnblocked, sumOp<label>()) << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::refinementHistory::apply
(
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    // Find common parent for all cells
    labelList cellToCluster;
    label nClusters = markCommonCells(cellToCluster);

    // Unblock all faces inbetween same cluster


    labelList clusterToProc(nClusters, -1);

    label nChanged = 0;

    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = mesh.faceOwner()[faceI];
        label nei = mesh.faceNeighbour()[faceI];

        label ownCluster = cellToCluster[own];
        label neiCluster = cellToCluster[nei];

        if (ownCluster != -1 && ownCluster == neiCluster)
        {
            if (clusterToProc[ownCluster] == -1)
            {
                clusterToProc[ownCluster] = decomposition[own];
            }

            if (decomposition[own] != clusterToProc[ownCluster])
            {
                decomposition[own] = clusterToProc[ownCluster];
                nChanged++;
            }
            if (decomposition[nei] != clusterToProc[ownCluster])
            {
                decomposition[nei] = clusterToProc[ownCluster];
                nChanged++;
            }
        }
    }

    if (refinementHistory::debug)
    {
        Info<< type() << " : changed decomposition on "
            << returnReduce(nChanged, sumOp<label>()) << " cells" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistory::refinementHistory(const IOobject& io)
:
    regIOobject(io),
    active_(false)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<refinementHistory>();

    readIOcontents();

    // When running in redistributePar + READ_IF_PRESENT it can happen
    // that some processors do have refinementHistory and some don't so
    // test for active has to be outside of above condition.
    active_ = returnReduceOr(visibleCells_.size());

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << " active:" << active_
            << endl;
    }
}


Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const List<splitCell8>& splitCells,
    const labelList& visibleCells,
    const bool active
)
:
    regIOobject(io),
    active_(active),
    splitCells_(splitCells),
    freeSplitCells_(),
    visibleCells_(visibleCells)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<refinementHistory>();

    readIOcontents();

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject or components :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << " active:" << active_
            << endl;
    }
}


Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const label nCells
)
:
    regIOobject(io),
    active_(false),
    freeSplitCells_()
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<refinementHistory>();

    if (!readIOcontents())
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell8());
        }
    }

    active_ = returnReduceOr(visibleCells_.size());


    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << " active:" << active_
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const label nCells,
    const bool active
)
:
    regIOobject(io),
    active_(active),
    freeSplitCells_()
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<refinementHistory>();

    if (!readIOcontents())
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label celli = 0; celli < nCells; celli++)
        {
            visibleCells_[celli] = celli;
            splitCells_.append(splitCell8());
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << " active:" << active_
            << endl;
    }
}


Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const refinementHistory& rh
)
:
    regIOobject(io),
    active_(rh.active_),
    splitCells_(rh.splitCells()),
    freeSplitCells_(rh.freeSplitCells()),
    visibleCells_(rh.visibleCells())
{
    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory : constructed initial"
            << " history." << endl;
    }
}


Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const UPtrList<const labelList>& cellMaps,
    const UPtrList<const refinementHistory>& refs
)
:
    regIOobject(io),
    active_(false)
{
    if (io.isAnyRead())
    {
        WarningInFunction
            << "Read option MUST_READ, READ_IF_PRESENT or READ_MODIFIED"
            << " suggests that a read constructor would be more appropriate."
            << endl;
    }

    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());


    // Determine offsets into splitCells
    labelList offsets(refs.size()+1);
    offsets[0] = 0;
    forAll(refs, refI)
    {
        const DynamicList<splitCell8>& subSplits = refs[refI].splitCells();
        offsets[refI+1] = offsets[refI]+subSplits.size();
    }

    // Construct merged splitCells
    splitCells_.setSize(offsets.last());
    forAll(refs, refI)
    {
        const DynamicList<splitCell8>& subSplits = refs[refI].splitCells();
        forAll(subSplits, i)
        {
            splitCell8& newSplit = splitCells_[offsets[refI]+i];

            // Copy
            newSplit = subSplits[i];

            // Offset indices
            if (newSplit.parent_ >= 0)
            {
                newSplit.parent_ += offsets[refI];
            }

            if (newSplit.addedCellsPtr_)
            {
                FixedList<label, 8>& splits = newSplit.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] += offsets[refI];
                    }
                }
            }
        }
    }


    // Construct merged visibleCells
    visibleCells_.setSize(mesh.nCells(), -1);
    forAll(refs, refI)
    {
        const labelList& cellMap = cellMaps[refI];
        const labelList& subVis = refs[refI].visibleCells();

        forAll(subVis, i)
        {
            label& newVis = visibleCells_[cellMap[i]];

            newVis = subVis[i];
            if (newVis >= 0)
            {
                newVis += offsets[refI];
            }
        }
    }


    // Is active if any of the refinementHistories is active (assumes active
    // flag parallel synchronised)
    active_ = false;
    forAll(refs, refI)
    {
        if (refs[refI].active())
        {
            active_ = true;
            break;
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from multiple refinementHistories :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


Foam::refinementHistory::refinementHistory(const IOobject& io, Istream& is)
:
    regIOobject(io),
    splitCells_(is),
    freeSplitCells_(),
    visibleCells_(is)
{
    active_ = returnReduceOr(visibleCells_.size());

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from Istream"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::refinementHistory> Foam::refinementHistory::clone
(
    const IOobject& io,
    // Per visible cell the processor it is going to
    const labelList& decomposition,
    // Per splitCell entry the processor it moves to
    const labelList& splitCellProc,
    // Per splitCell entry the number of live cells that move to that processor
    const labelList& splitCellNum,

    const label procI,

    // From old to new splitCells
    labelList& oldToNewSplit
) const
{
    oldToNewSplit.setSize(splitCells_.size());
    oldToNewSplit = -1;

    // Compacted splitCells
    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // Loop over all entries. Note: could recurse like countProc so only
    // visit used entries but is probably not worth it.

    forAll(splitCells_, index)
    {
        if (splitCellProc[index] == procI && splitCellNum[index] == 8)
        {
            // Entry moves in its whole to procI
            oldToNewSplit[index] = newSplitCells.size();
            newSplitCells.append(splitCells_[index]);
        }
    }

    // Add live cells that are subsetted.
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0 && decomposition[cellI] == procI)
        {
            label parent = splitCells_[index].parent_;

            // Create new splitCell with parent
            oldToNewSplit[index] = newSplitCells.size();
            newSplitCells.append(splitCell8(parent));
        }
    }

    //forAll(oldToNewSplit, index)
    //{
    //    Pout<< "old:" << index << " new:" << oldToNewSplit[index]
    //        << endl;
    //}

    newSplitCells.shrink();

    // Renumber contents of newSplitCells
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNewSplit[split.parent_];
        }
        if (split.addedCellsPtr_)
        {
            FixedList<label, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    splits[i] = oldToNewSplit[splits[i]];
                }
            }
        }
    }


    // Count number of cells
    label nSub = 0;
    forAll(decomposition, cellI)
    {
        if (decomposition[cellI] == procI)
        {
            nSub++;
        }
    }

    labelList newVisibleCells(nSub);
    nSub = 0;

    forAll(visibleCells_, cellI)
    {
        if (decomposition[cellI] == procI)
        {
            label index = visibleCells_[cellI];
            if (index >= 0)
            {
                index = oldToNewSplit[index];
            }
            newVisibleCells[nSub++] = index;
        }
    }

    return autoPtr<refinementHistory>
    (
        new refinementHistory
        (
            io,
            newSplitCells,
            newVisibleCells,
            active_
        )
    );
}


Foam::autoPtr<Foam::refinementHistory> Foam::refinementHistory::clone
(
    const IOobject& io,
    const labelList& cellMap
) const
{
    if (active_)
    {
        // Mark selected cells with '1'
        labelList decomposition(visibleCells_.size(), Zero);
        forAll(cellMap, i)
        {
            decomposition[cellMap[i]] = 1;
        }


        // Per splitCell entry the processor it moves to
        labelList splitCellProc(splitCells_.size(), -1);
        // Per splitCell entry the number of live cells that move to that
        // processor
        labelList splitCellNum(splitCells_.size(), Zero);

        forAll(visibleCells_, cellI)
        {
            label index = visibleCells_[cellI];

            if (index >= 0)
            {
                countProc
                (
                    splitCells_[index].parent_,
                    decomposition[cellI],
                    splitCellProc,
                    splitCellNum
                );
            }
        }

        labelList oldToNewSplit;
        return clone
        (
            io,
            decomposition,
            splitCellProc,
            splitCellNum,
            1,      //procI,
            oldToNewSplit
        );
    }
    else
    {
        return autoPtr<refinementHistory>
        (
            new refinementHistory
            (
                io,
                List<splitCell8>(),
                labelList(),
                false
            )
        );
    }
}


void Foam::refinementHistory::resize(const label size)
{
    label oldSize = visibleCells_.size();

    if (debug)
    {
        Pout<< "refinementHistory::resize from " << oldSize << " to " << size
            << " cells" << endl;
    }

    visibleCells_.setSize(size);

    // Set additional elements to -1.
    for (label i = oldSize; i < visibleCells_.size(); i++)
    {
        visibleCells_[i] = -1;
    }
}


void Foam::refinementHistory::updateMesh(const mapPolyMesh& map)
{
    if (active())
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        // Note that only the live cells need to be renumbered.

        labelList newVisibleCells(map.cellMap().size(), -1);

        forAll(visibleCells_, celli)
        {
            if (visibleCells_[celli] != -1)
            {
                label index = visibleCells_[celli];

                // Check not already set
                if (splitCells_[index].addedCellsPtr_)
                {
                    FatalErrorInFunction
                        << "Problem" << abort(FatalError);
                }

                label newCelli = reverseCellMap[celli];

                if (newCelli >= 0)
                {
                    newVisibleCells[newCelli] = index;
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementHistory::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    if (active())
    {
        labelList newVisibleCells(cellMap.size(), -1);

        forAll(newVisibleCells, celli)
        {
            label oldCelli = cellMap[celli];

            label index = visibleCells_[oldCelli];

            // Check that cell is live (so its parent has no refinement)
            if (index >= 0 && splitCells_[index].addedCellsPtr_)
            {
                FatalErrorInFunction
                    << "Problem" << abort(FatalError);
            }

            newVisibleCells[celli] = index;
        }

        if (debug)
        {
            Pout<< "refinementHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementHistory::countProc
(
    const label index,
    const label newProcNo,
    labelList& splitCellProc,
    labelList& splitCellNum
) const
{
    if (splitCellProc[index] != newProcNo)
    {
        // Different destination processor from other cells using this
        // parent. Reset count.
        splitCellProc[index] = newProcNo;
        splitCellNum[index] = 1;
    }
    else
    {
        splitCellNum[index]++;

        // Increment parent if whole splitCell moves to same processor
        if (splitCellNum[index] == 8)
        {
            if (debug)
            {
                Pout<< "Moving " << splitCellNum[index]
                    << " cells originating from cell " << index
                    << " from processor " << Pstream::myProcNo()
                    << " to processor " << splitCellProc[index]
                    << endl;
            }

            label parent = splitCells_[index].parent_;

            if (parent >= 0)
            {
                countProc(parent, newProcNo, splitCellProc, splitCellNum);
            }
        }
    }
}


void Foam::refinementHistory::distribute(const mapDistributePolyMesh& map)
{
    if (!active())
    {
        FatalErrorInFunction
            << "Calling distribute on inactive history" << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        return;
    }

    // Remove unreferenced history.
    compact();

    //Pout<< nl << "--BEFORE:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;


    // Distribution is only partially functional.
    // If all 8 cells resulting from a single parent are sent across in one
    // go it will also send across that part of the refinement history.
    // If however e.g. first 1 and then the other 7 are sent across the
    // history will not be reconstructed.

    // Determine clusters. This is per every entry in splitCells_ (that is
    // a parent of some refinement) a label giving the processor it goes to
    // if all its children are going to the same processor.

    // Per visible cell the processor it goes to.
    labelList destination(visibleCells_.size());

    const labelListList& subCellMap = map.cellMap().subMap();

    forAll(subCellMap, proci)
    {
        const labelList& newToOld = subCellMap[proci];

        forAll(newToOld, i)
        {
            label oldCelli = newToOld[i];

            destination[oldCelli] = proci;
        }
    }

    // Per splitCell entry the processor it moves to
    labelList splitCellProc(splitCells_.size(), -1);
    // Per splitCell entry the number of live cells that move to that processor
    labelList splitCellNum(splitCells_.size(), Zero);

    forAll(visibleCells_, celli)
    {
        label index = visibleCells_[celli];

        if (index >= 0)
        {
            countProc
            (
                splitCells_[index].parent_,
                destination[celli],
                splitCellProc,
                splitCellNum
            );
        }
    }

    //Pout<< "refinementHistory::distribute :"
    //    << " splitCellProc:" << splitCellProc << endl;
    //
    //Pout<< "refinementHistory::distribute :"
    //    << " splitCellNum:" << splitCellNum << endl;


    // Create subsetted refinement tree consisting of all parents that
    // move in their whole to other processor.

    PstreamBuffers pBufs;

    for (const int proci : Pstream::allProcs())
    {
        //Pout<< "-- Subetting for processor " << proci << endl;

        // From uncompacted to compacted splitCells.
        labelList oldToNew(splitCells_.size(), -1);

        // Compacted splitCells. Similar to subset routine below.
        DynamicList<splitCell8> newSplitCells(splitCells_.size());

        // Loop over all entries. Note: could recurse like countProc so only
        // visit used entries but is probably not worth it.

        forAll(splitCells_, index)
        {
            if (splitCellProc[index] == proci && splitCellNum[index] == 8)
            {
                // Entry moves in its whole to proci
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCells_[index]);
            }
        }

        // Add live cells that are subsetted.
        forAll(visibleCells_, celli)
        {
            label index = visibleCells_[celli];

            if (index >= 0 && destination[celli] == proci)
            {
                label parent = splitCells_[index].parent_;

                // Create new splitCell with parent
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCell8(parent));
            }
        }

        //forAll(oldToNew, index)
        //{
        //    Pout<< "old:" << index << " new:" << oldToNew[index]
        //        << endl;
        //}

        newSplitCells.shrink();

        // Renumber contents of newSplitCells
        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ = oldToNew[split.parent_];
            }
            if (split.addedCellsPtr_)
            {
                FixedList<label, 8>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] = oldToNew[splits[i]];
                    }
                }
            }
        }


        const labelList& subMap = subCellMap[proci];

        // New visible cells.
        labelList newVisibleCells(subMap.size(), -1);

        forAll(subMap, newCelli)
        {
            label oldCelli = subMap[newCelli];

            label oldIndex = visibleCells_[oldCelli];

            if (oldIndex >= 0)
            {
                newVisibleCells[newCelli] = oldToNew[oldIndex];
            }
        }

        //Pout<< nl << "--Subset for domain:" << proci << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // Send to neighbours
        UOPstream toNbr(proci, pBufs);
        toNbr << newSplitCells << newVisibleCells;
    }


    // Wait for finish
    pBufs.finishedSends();


    // Receive from neighbours and merge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Remove all entries. Leave storage intact.
    splitCells_.clear();

    const polyMesh& mesh = dynamic_cast<const polyMesh&>(db());

    visibleCells_.setSize(mesh.nCells());
    visibleCells_ = -1;

    for (const int proci : Pstream::allProcs())
    {
        UIPstream fromNbr(proci, pBufs);
        List<splitCell8> newSplitCells(fromNbr);
        labelList newVisibleCells(fromNbr);

        //Pout<< nl << "--Received from domain:" << proci << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // newSplitCells contain indices only into newSplitCells so
        // renumbering can be done here.
        label offset = splitCells_.size();

        //Pout<< "**Renumbering data from proc " << proci << " with offset "
        //    << offset << endl;

        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ += offset;
            }
            if (split.addedCellsPtr_)
            {
                FixedList<label, 8>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] += offset;
                    }
                }
            }

            splitCells_.append(split);
        }


        // Combine visibleCell.
        const labelList& constructMap = map.cellMap().constructMap()[proci];

        forAll(newVisibleCells, i)
        {
            if (newVisibleCells[i] >= 0)
            {
                visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
            }
        }
    }
    splitCells_.shrink();

    //Pout<< nl << "--AFTER:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;
}


void Foam::refinementHistory::compact()
{
    if (debug)
    {
        Pout<< "refinementHistory::compact() Entering with:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;

        // Check all free splitCells are marked as such
        forAll(freeSplitCells_, i)
        {
            label index = freeSplitCells_[i];

            if (splitCells_[index].parent_ != -2)
            {
                FatalErrorInFunction
                    << "Problem index:" << index
                    << abort(FatalError);
            }
        }

        // Check none of the visible cells are marked as free
        forAll(visibleCells_, celli)
        {
            if
            (
                visibleCells_[celli] >= 0
             && splitCells_[visibleCells_[celli]].parent_ == -2
            )
            {
                FatalErrorInFunction
                    << "Problem : visible cell:" << celli
                    << " is marked as being free." << abort(FatalError);
            }
        }
    }

    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // From uncompacted to compacted splitCells.
    labelList oldToNew(splitCells_.size(), -1);

    // Mark all used splitCell entries. These are either indexed by visibleCells
    // or indexed from other splitCell entries.

    // Mark from visibleCells
    forAll(visibleCells_, celli)
    {
        label index = visibleCells_[celli];

        if (index >= 0)
        {
            // Make sure we only mark visible indices if they either have a
            // parent or subsplits.
            if
            (
                splitCells_[index].parent_ != -1
             || splitCells_[index].addedCellsPtr_
            )
            {
                markSplit(index, oldToNew, newSplitCells);
            }
        }
    }

    // Mark from splitCells
    forAll(splitCells_, index)
    {
        if (splitCells_[index].parent_ == -2)
        {
            // freed cell.
        }
        else if
        (
            splitCells_[index].parent_ == -1
         && !splitCells_[index].addedCellsPtr_
        )
        {
            // recombined cell. No need to keep since no parent and no subsplits
            // Note that gets marked if reachable from other index!
        }
        else
        {
            // Is used element.
            markSplit(index, oldToNew, newSplitCells);
        }
    }


    // Now oldToNew is fully complete and compacted elements are in
    // newSplitCells.
    // Renumber contents of newSplitCells and visibleCells.
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNew[split.parent_];
        }
        if (split.addedCellsPtr_)
        {
            FixedList<label, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    splits[i] = oldToNew[splits[i]];
                }
            }
        }
    }


    if (debug)
    {
        Pout<< "refinementHistory::compact : compacted splitCells from "
            << splitCells_.size() << " to " << newSplitCells.size() << endl;
    }

    splitCells_.transfer(newSplitCells);
    freeSplitCells_.clearStorage();


    if (debug)
    {
        Pout<< "refinementHistory::compact() NOW:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " newSplitCells:" << newSplitCells.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
    }


    // Adapt indices in visibleCells_
    forAll(visibleCells_, celli)
    {
        label index = visibleCells_[celli];

        if (index >= 0)
        {
            // Note that oldToNew can be -1 so it resets newVisibleCells.
            visibleCells_[celli] = oldToNew[index];
        }
        else
        {
            // Keep -1 value.
        }
    }
}


void Foam::refinementHistory::writeDebug() const
{
    writeDebug(visibleCells_, splitCells_);
}


void Foam::refinementHistory::storeSplit
(
    const label celli,
    const labelList& addedCells
)
{
    label parentIndex = -1;

    if (visibleCells_[celli] != -1)
    {
        // Was already live. The current live cell becomes the
        // parent of the cells split off from it.

        parentIndex = visibleCells_[celli];

        // It is no longer live (note that actually celli gets alive
        // again below since is addedCells[0])
        visibleCells_[celli] = -1;
    }
    else
    {
        // Create 0th level. -1 parent to denote this.
        parentIndex = allocateSplitCell(-1, -1);
    }

    // Create live entries for added cells that point to the
    // cell they were created from (parentIndex)
    forAll(addedCells, i)
    {
        label addedCelli = addedCells[i];

        // Create entries for the split off cells. All of them
        // are visible.
        visibleCells_[addedCelli] = allocateSplitCell(parentIndex, i);
    }
}


void Foam::refinementHistory::combineCells
(
    const label masterCelli,
    const labelList& combinedCells
)
{
    // Save the parent structure
    label parentIndex = splitCells_[visibleCells_[masterCelli]].parent_;

    // Remove the information for the combined cells
    forAll(combinedCells, i)
    {
        label celli = combinedCells[i];

        freeSplitCell(visibleCells_[celli]);
        visibleCells_[celli] = -1;
    }

    splitCell8& parentSplit = splitCells_[parentIndex];
    parentSplit.addedCellsPtr_.reset(nullptr);
    visibleCells_[masterCelli] = parentIndex;
}


bool Foam::refinementHistory::read()
{
    bool ok = readData(readStream(typeName));
    close();

    active_ = returnReduceOr(visibleCells_.size());

    return ok;
}


bool Foam::refinementHistory::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::refinementHistory::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


void Foam::refinementHistory::removeFiles(const polyMesh& mesh)
{
    IOobject io
    (
        "dummy",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh
    );
    fileName setsDir(io.path());

    if (topoSet::debug) DebugVar(setsDir);

    if (exists(setsDir/typeName))
    {
        rm(setsDir/typeName);
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistory& rh)
{
    rh.freeSplitCells_.clearStorage();

    is >> rh.splitCells_ >> rh.visibleCells_;

    // Check indices.
    rh.checkIndices();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementHistory& rh)
{
    const_cast<refinementHistory&>(rh).compact();

    return os   << "// splitCells" << nl
                << rh.splitCells_ << nl
                << "// visibleCells" << nl
                << rh.visibleCells_;
}


// ************************************************************************* //
