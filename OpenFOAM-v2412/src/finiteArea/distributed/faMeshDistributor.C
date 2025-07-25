/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "faMeshDistributor.H"
#include "BitOps.H"
#include "fileOperation.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::faMeshDistributor::verbose_ = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMeshDistributor::createPatchMaps() const
{
    const faBoundaryMesh& oldPatches = srcMesh_.boundary();
    const faBoundaryMesh& newPatches = tgtMesh_.boundary();

    patchEdgeMaps_.clear();
    patchEdgeMaps_.resize(oldPatches.size());

    // area: edgeMap (volume: faceMap)
    const auto& faEdgeMap = distMap_.faceMap();

    // The logical edge ranges per patch [target mesh]
    List<labelRange> ranges = newPatches.patchRanges();

    forAll(oldPatches, patchi)
    {
        if (!isA<processorFaPatch>(oldPatches[patchi]))
        {
            // Map non-processor only

            // Copy full map
            patchEdgeMaps_.set
            (
                patchi,
                new mapDistributeBase(faEdgeMap)
            );

            // Retain patch elements
            patchEdgeMaps_[patchi].compactRemoteData
            (
                bitSet(ranges[patchi]),
                UPstream::msgType(),
                true  // Also renumber/resize the compact maps
            );
        }
    }
}


void Foam::faMeshDistributor::createInternalEdgeMap() const
{
    // area: edgeMap (volume: faceMap)
    const auto& faEdgeMap = distMap_.faceMap();

    // Copy full map
    internalEdgeMapPtr_ = std::make_unique<mapDistributeBase>(faEdgeMap);

    // Retain internal edges
    internalEdgeMapPtr_->compactRemoteData
    (
        bitSet(tgtMesh_.nInternalEdges(), true),
        UPstream::msgType(),
        true  // Also renumber/resize the compact maps
    );
}


void Foam::faMeshDistributor::checkAddressing() const
{
    #ifdef FULLDEBUG
    {
        Pout<< "Create from nFaces:" << srcMesh_.faceLabels().size()
            << " to:" << tgtMesh_.faceLabels().size() << endl;

        // Check face centres
        {
            vectorField oldFaceCentres(srcMesh_.areaCentres());
            vectorField newFaceCentres(tgtMesh_.areaCentres());

            // volume: cells, area: faces
            distMap_.distributeCellData(oldFaceCentres);

            vectorField diff(newFaceCentres - oldFaceCentres);

            if (!diff.empty() && !diff.uniform())
            {
                forAll(oldFaceCentres, facei)
                {
                    if (oldFaceCentres[facei] != newFaceCentres[facei])
                    {
                        Pout<< "face: " << facei
                            << ' ' << oldFaceCentres[facei]
                            << " vs "  << newFaceCentres[facei]
                            << endl;
                    }
                }
            }
        }

        // Check edge centres
        {
            vectorField oldEdgeCentres
            (
                faMeshTools::flattenEdgeField(srcMesh_.edgeCentres())
            );
            vectorField newEdgeCentres
            (
                faMeshTools::flattenEdgeField(tgtMesh_.edgeCentres())
            );

            Pout<< "distributed edges: " << oldEdgeCentres.size() << " from "
                << srcMesh_.nEdges() << " to " << tgtMesh_.nEdges() << endl;

            // volume: faces, area: edges
            distMap_.distributeFaceData(oldEdgeCentres);

            vectorField diff(newEdgeCentres - oldEdgeCentres);

            if (!diff.empty() && !diff.uniform())
            {
                forAll(oldEdgeCentres, edgei)
                {
                    if (oldEdgeCentres[edgei] != newEdgeCentres[edgei])
                    {
                        Pout<< "edge: " << edgei
                            << ' ' << oldEdgeCentres[edgei]
                            << " vs "  << newEdgeCentres[edgei]
                            << endl;
                    }
                }
            }

            Info<< "Patch edge maps" << endl;
            forAll(patchEdgeMaps_, patchi)
            {
                if (patchEdgeMaps_.set(patchi))
                {
                    Pout<< "patch " << patchi << " : "
                        << patchEdgeMaps_[patchi].info() << endl;
                }
            }

            Info<< nl << "Detailed patch maps" << endl;

            forAll(patchEdgeMaps_, patchi)
            {
                if (patchEdgeMaps_.set(patchi))
                {
                    Info<< "patch " << patchi << " : "
                        << patchEdgeMaps_[patchi] << endl;
                }
            }
        }
    }
    #endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshDistributor::faMeshDistributor
(
    const faMesh& srcMesh,
    const faMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    const bool isWriteProc
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap),
    dummyHandler_(fileOperation::null()),
    writeHandler_(dummyHandler_),
    isWriteProc_(isWriteProc)
{
    checkAddressing();
}


Foam::faMeshDistributor::faMeshDistributor
(
    const faMesh& srcMesh,
    const faMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    refPtr<fileOperation>& writeHandler
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap),
    dummyHandler_(nullptr),
    writeHandler_(writeHandler),
    isWriteProc_(Switch::INVALID)
{
    checkAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faMeshDistributor::distributeAllFields
(
    const IOobjectList& objects,
    const wordRes& selected
) const
{
    label nTotal = 0;

    nTotal += distributeAreaFields<scalar>(objects, selected);
    nTotal += distributeAreaFields<vector>(objects, selected);
    nTotal += distributeAreaFields<symmTensor>(objects, selected);
    nTotal += distributeAreaFields<sphericalTensor>(objects, selected);
    nTotal += distributeAreaFields<tensor>(objects, selected);

    nTotal += distributeEdgeFields<scalar>(objects, selected);
    nTotal += distributeEdgeFields<vector>(objects, selected);
    nTotal += distributeEdgeFields<symmTensor>(objects, selected);
    nTotal += distributeEdgeFields<sphericalTensor>(objects, selected);
    nTotal += distributeEdgeFields<tensor>(objects, selected);

    return nTotal;
}


// ************************************************************************* //
