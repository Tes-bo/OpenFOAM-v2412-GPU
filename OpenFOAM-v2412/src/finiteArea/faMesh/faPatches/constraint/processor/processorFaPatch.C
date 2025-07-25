/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "processorFaPatch.H"
#include "processorPolyPatch.H"  // For newName()
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "faMesh.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorFaPatch, 0);
    addToRunTimeSelectionTable(faPatch, processorFaPatch, dictionary);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorFaPatch::~processorFaPatch()
{
    neighbPointsPtr_.clear();
    nonGlobalPatchPointsPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFaPatch::processorFaPatch
(
    const word& name,
    const labelUList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label nbrPolyPatchi,
    const label myProcNo,
    const label neighbProcNo,
    const word& patchType
)
:
    coupledFaPatch(name, edgeLabels, index, bm, nbrPolyPatchi, patchType),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbEdgeCentres_(),
    neighbEdgeLengths_(),
    neighbEdgeFaceCentres_(),
    neighbPointsPtr_(nullptr),
    nonGlobalPatchPointsPtr_(nullptr)
{}


Foam::processorFaPatch::processorFaPatch
(
    const labelUList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label nbrPolyPatchi,
    const label myProcNo,
    const label neighbProcNo,
    const word& patchType
)
:
    processorFaPatch
    (
        processorPolyPatch::newName(myProcNo, neighbProcNo),
        edgeLabels,
        index,
        bm,
        nbrPolyPatchi,
        myProcNo,
        neighbProcNo,
        patchType
    )
{}


Foam::processorFaPatch::processorFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm,
    const word& patchType
)
:
    coupledFaPatch(name, dict, index, bm, patchType),
    myProcNo_(dict.get<label>("myProcNo")),
    neighbProcNo_(dict.get<label>("neighbProcNo")),
    neighbEdgeCentres_(),
    neighbEdgeLengths_(),
    neighbEdgeFaceCentres_(),
    neighbPointsPtr_(nullptr),
    nonGlobalPatchPointsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::processorFaPatch::comm() const
{
    return boundaryMesh().mesh().comm();
}


void Foam::processorFaPatch::makeNonGlobalPatchPoints() const
{
    // If it is not running parallel or there are no global points
    // create a 1->1 map

    // Can not use faGlobalMeshData at this point yet
    // - use polyMesh globalData instead

    const auto& pMeshGlobalData = boundaryMesh().mesh().mesh().globalData();

    if (!Pstream::parRun() || !pMeshGlobalData.nGlobalPoints())
    {
        // 1 -> 1 mapping
        nonGlobalPatchPointsPtr_.reset(new labelList(identity(nPoints())));
    }
    else
    {
        nonGlobalPatchPointsPtr_.reset(new labelList(nPoints()));
        labelList& ngpp = *nonGlobalPatchPointsPtr_;

        // Get reference to shared points
        const labelList& sharedPoints = pMeshGlobalData.sharedPointLabels();

        const labelList& faMeshPatchPoints = pointLabels();

        const labelList& meshPoints =
            boundaryMesh().mesh().patch().meshPoints();

        label nNonShared = 0;

        forAll(faMeshPatchPoints, pointi)
        {
            const label mpi = meshPoints[faMeshPatchPoints[pointi]];
            if (!sharedPoints.found(mpi))
            {
                ngpp[nNonShared] = pointi;
                ++nNonShared;
            }
        }

        ngpp.resize(nNonShared);
    }
}


void Foam::processorFaPatch::initGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        if (neighbProcNo() >= pBufs.nProcs())
        {
            FatalErrorInFunction
                << "On patch " << name()
                << " trying to access out of range neighbour processor "
                << neighbProcNo() << ". This can happen if" << nl
                << "    trying to run on an incorrect number of processors"
                << exit(FatalError);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << edgeCentres()
            << edgeLengths()
            << edgeFaceCentres();
    }
}


void Foam::processorFaPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> neighbEdgeCentres_
                >> neighbEdgeLengths_
                >> neighbEdgeFaceCentres_;
        }

        #ifdef FULLDEBUG
        const scalarField& magEl = magEdgeLengths();

        forAll(magEl, edgei)
        {
            scalar nmagEl = mag(neighbEdgeLengths_[edgei]);
            scalar avEl = (magEl[edgei] + nmagEl)/2.0;

            if (mag(magEl[edgei] - nmagEl)/avEl > 1e-6)
            {
                WarningInFunction
                    << "edge " << edgei
                    << " length does not match neighbour by "
                    << 100*mag(magEl[edgei] - nmagEl)/avEl
                    << "% -- possible edge ordering problem" << nl;
            }
        }
        #endif

        calcTransformTensors
        (
            edgeCentres(),
            neighbEdgeCentres_,
            edgeNormals(),
            neighbEdgeNormals()
        );
    }
}


void Foam::processorFaPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    faPatch::movePoints(pBufs, p);
    initGeometry(pBufs);
}


void Foam::processorFaPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    processorFaPatch::calcGeometry(pBufs);
}


void Foam::processorFaPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    faPatch::initUpdateMesh(pBufs);

    neighbPointsPtr_.clear();

    if (Pstream::parRun())
    {
        if (neighbProcNo() >= pBufs.nProcs())
        {
            FatalErrorInFunction
                << "On patch " << name()
                << " trying to access out of range neighbour processor "
                << neighbProcNo() << ". This can happen if" << nl
                << "    trying to run on an incorrect number of processors"
                << exit(FatalError);
        }

        // Express all points as patch edge and index in edge.
        labelList patchEdge(nPoints());
        labelList indexInEdge(nPoints());

        const edgeList::subList patchEdges =
            patchSlice(boundaryMesh().mesh().edges());

        const labelListList& ptEdges = pointEdges();

        for (label patchPointI = 0; patchPointI < nPoints(); ++patchPointI)
        {
            label edgeI = ptEdges[patchPointI][0];

            patchEdge[patchPointI] = edgeI;

            const edge& e = patchEdges[edgeI];

            indexInEdge[patchPointI] = e.find(pointLabels()[patchPointI]);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << patchEdge
            << indexInEdge;
    }
}


void Foam::processorFaPatch::updateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    faPatch::updateMesh(pBufs);

    neighbPointsPtr_.clear();

    if (Pstream::parRun())
    {
        labelList nbrPatchEdge;
        labelList nbrIndexInEdge;

        {
            // Note cannot predict exact size since edgeList not (yet) sent as
            // binary entity but as List of edges.

            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrPatchEdge
                >> nbrIndexInEdge;
        }

        if (nbrPatchEdge.size() == nPoints())
        {
            // Convert neighbour edges and indices into face back into
            // my edges and points.
            neighbPointsPtr_.reset(new labelList(nPoints()));
            labelList& neighbPoints = *neighbPointsPtr_;

            const edgeList::subList patchEdges =
                patchSlice(boundaryMesh().mesh().edges());

            forAll(nbrPatchEdge, nbrPointI)
            {
                // Find edge and index in edge on this side.
                const edge& e = patchEdges[nbrPatchEdge[nbrPointI]];

                const label index = 1 - nbrIndexInEdge[nbrPointI];

                const label patchPointI = pointLabels().find(e[index]);

                neighbPoints[patchPointI] = nbrPointI;
            }
        }
        else
        {
            // Differing number of points. Probably patch includes
            // part of a cyclic.
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::processorFaPatch::neighbEdgeNormals() const
{
    auto tresult = tmp<vectorField>::New(neighbEdgeLengths_);
    tresult.ref().normalise();
    return tresult;
}


const Foam::labelList& Foam::processorFaPatch::neighbPoints() const
{
    if (!neighbPointsPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points  on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }

    return *neighbPointsPtr_;
}


void Foam::processorFaPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        const vectorField& skew = skewCorrectionVectors();
        const scalarField& PN = lPN();

        const scalarField lEN
        (
            mag(neighbEdgeFaceCentres() - edgeCentres() + skew)
        );

        forAll(w, i)
        {
            if (mag(PN[i]) > SMALL)
            {
                w[i] = lEN[i]/PN[i];
            }
        }
    }
    else
    {
        w = scalar(1);
    }
}


void Foam::processorFaPatch::makeLPN(scalarField& lPN) const
{
    const vectorField& skew = skewCorrectionVectors();

    if (Pstream::parRun())
    {
        lPN =
            mag(edgeCentres() - skew - edgeFaceCentres())         // lPE
          + mag(neighbEdgeFaceCentres() - edgeCentres() + skew);  // lEN

        for (scalar& lpn: lPN)
        {
            if (mag(lpn) < SMALL)
            {
                lpn = SMALL;
            }
        }
    }
    else
    {
        lPN =
            mag(edgeCentres() - skew - edgeFaceCentres())         // lPE
          + mag(edgeFaceCentres() - edgeCentres() + skew);        // lEN
    }
}


void Foam::processorFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (Pstream::parRun())
    {
        const edgeList& edges = boundaryMesh().mesh().edges();
        const pointField& points = boundaryMesh().mesh().points();
        const vectorField& lengths = edgeLengths();
        const scalarField& PN = lPN();

        tmp<vectorField> tdelta = processorFaPatch::delta();
        vectorField& unitDelta = tdelta.ref();

        forAll(dc, i)
        {
            vector edgeNormal =
                normalised(lengths[i] ^ edges[i + start()].vec(points));

            unitDelta[i].removeCollinear(edgeNormal);
            unitDelta[i].normalise();

            edgeNormal = normalised(lengths[i]);

            const scalar alpha = PN[i]*(edgeNormal & unitDelta[i]);
            if (mag(alpha) > SMALL)
            {
                dc[i] = scalar(1)/max(alpha, 0.05*PN[i]);
            }
        }
    }
    else
    {
        dc = scalar(1)
            /((edgeNormals() & coupledFaPatch::delta()) + VSMALL);
    }
}


void Foam::processorFaPatch::makeCorrectionVectors(vectorField& cv) const
{
    if (Pstream::parRun())
    {
        const edgeList& edges = boundaryMesh().mesh().edges();
        const pointField& points = boundaryMesh().mesh().points();
        const vectorField& lengths = edgeLengths();

        tmp<vectorField> tdelta = processorFaPatch::delta();
        vectorField& unitDelta = tdelta.ref();

        forAll(cv, i)
        {
            vector edgeNormal =
                normalised(lengths[i] ^ edges[i + start()].vec(points));

            unitDelta[i].removeCollinear(edgeNormal);
            unitDelta.normalise();

            edgeNormal = normalised(lengths[i]);

            const scalar alpha = unitDelta[i] & edgeNormal;
            scalar dc = SMALL;
            if (mag(alpha) > SMALL)
            {
                dc = scalar(1)/alpha;
            }

            cv[i] = edgeNormal - dc*unitDelta[i];
        }
    }
    else
    {
        cv = Zero;
    }
}


Foam::tmp<Foam::vectorField> Foam::processorFaPatch::delta() const
{
    if (Pstream::parRun())
    {
        // Do the transformation if necessary
        if (parallel())
        {
            return
                coupledFaPatch::delta()
              - (
                    neighbEdgeCentres()
                  - neighbEdgeFaceCentres()
                );
        }
        else
        {
            return
                coupledFaPatch::delta()
              - transform
                (
                    forwardT(),
                    (
                        neighbEdgeCentres()
                      - neighbEdgeFaceCentres()
                    )
                );
        }
    }
    else
    {
        return coupledFaPatch::delta();
    }
}


const Foam::labelList& Foam::processorFaPatch::nonGlobalPatchPoints() const
{
    if (!nonGlobalPatchPointsPtr_)
    {
        makeNonGlobalPatchPoints();
    }

    return *nonGlobalPatchPointsPtr_;
}


Foam::tmp<Foam::labelField> Foam::processorFaPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::processorFaPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& edgeFaces
) const
{
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, edgeFaces, tpfld.ref());
    return tpfld;
}


void Foam::processorFaPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    send(commsType, interfaceData);
}


Foam::tmp<Foam::labelField> Foam::processorFaPatch::transfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorFaPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


Foam::tmp<Foam::labelField> Foam::processorFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorFaPatch::write(Ostream& os) const
{
    faPatch::write(os);
    os.writeEntry("myProcNo",  myProcNo_);
    os.writeEntry("neighbProcNo", neighbProcNo_);
}


// ************************************************************************* //
