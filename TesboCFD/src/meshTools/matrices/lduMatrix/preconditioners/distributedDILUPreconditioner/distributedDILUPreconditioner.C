/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "distributedDILUPreconditioner.H"
#include "processorLduInterface.H"
#include "cyclicLduInterface.H"
#include "cyclicAMILduInterface.H"
#include "processorColour.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedDILUPreconditioner, 0);

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<distributedDILUPreconditioner>
        adddistributedDILUPreconditionerAsymMatrixConstructorToTable_;

    const lduMesh* distributedDILUPreconditioner::meshPtr_ = nullptr;

    autoPtr<labelList> distributedDILUPreconditioner::procColoursPtr_(nullptr);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributedDILUPreconditioner::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& coupleCoeffs,
    const labelList& selectedInterfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();
    const auto& interfaces = solver_.interfaces();

    const label startOfRequests = Pstream::nRequests();

    // From lduMatrix::initMatrixInterfaces

    // Avoid any conflicts with inter-processor comms
    const int oldTag = UPstream::incrMsgType(321);

    for (const label interfacei : selectedInterfaces)
    {
        interfaces[interfacei].initInterfaceMatrixUpdate
        (
            result,
            add,
            lduAddr,
            interfacei,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            UPstream::commsTypes::nonBlocking
        );
    }

    UPstream::waitRequests(startOfRequests);

    // Consume
    for (const label interfacei : selectedInterfaces)
    {
        interfaces[interfacei].updateInterfaceMatrix
        (
            result,
            add,
            lduAddr,
            interfacei,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            UPstream::commsTypes::nonBlocking
        );
    }

    UPstream::msgType(oldTag);  // Restore tag
}


void Foam::distributedDILUPreconditioner::sendGlobal
(
    const labelList& selectedInterfaces,
    solveScalarField& psi,
    const label colouri
) const
{
    const auto& interfaces = solver_.interfaces();

    if (selectedInterfaces.size())
    {
        // Save old data
        FieldField<Field, scalar> one(interfaces.size());
        FieldField<Field, solveScalar> old(interfaces.size());
        for (const label inti : selectedInterfaces)
        {
            const auto& intf = interfaces[inti].interface();
            const auto& fc = intf.faceCells();
            one.set(inti, new scalarField(fc.size(), 1.0));
            old.set(inti, new solveScalarField(psi, intf.faceCells()));
        }
        updateMatrixInterfaces
        (
            false,              // add to psi
            one,
            selectedInterfaces,
            psi,                // send data
            psi,                // result
            0                   // cmpt
        );

        if (!colourBufs_.set(colouri))
        {
            colourBufs_.set
            (
                colouri,
                new FieldField<Field, solveScalar>
                (
                    interfaces.size()
                )
            );
        }
        auto& colourBuf = colourBufs_[colouri];
        colourBuf.setSize(interfaces.size());
        for (const label inti : selectedInterfaces)
        {
            const auto& intf = interfaces[inti].interface();
            const auto& fc = intf.faceCells();
            if (!colourBuf.set(inti))
            {
                colourBuf.set(inti, new Field<solveScalar>(fc.size()));
            }
            auto& cb = colourBuf[inti];
            //cb.resize_nocopy(fc.size());
            auto& oldValues = old[inti];

            forAll(cb, face)
            {
                const label cell = fc[face];
                // Store change in boundary value
                cb[face] = psi[cell]-oldValues[face];
                // Restore old value
                std::swap(psi[cell], oldValues[face]);
            }
        }
    }
}


void Foam::distributedDILUPreconditioner::receive
(
    const labelList& selectedInterfaces,
    DynamicList<UPstream::Request>& requests
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();

    // Start reads
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto* ppp = isA<processorLduInterface>(intf);

        auto& recvBuf = recvBufs_[inti];
        recvBuf.resize_nocopy(interfaceBouCoeffs[inti].size());

        requests.push_back(UPstream::Request());
        UIPstream::read
        (
            requests.back(),
            ppp->neighbProcNo(),
            recvBuf.data_bytes(),
            recvBuf.size_bytes(),
            ppp->tag()+70,          // random offset
            ppp->comm()
        );
    }
}


void Foam::distributedDILUPreconditioner::send
(
    const labelList& selectedInterfaces,
    const solveScalarField& psiInternal,
    DynamicList<UPstream::Request>& requests
) const
{
    const auto& interfaces = solver_.interfaces();

    // Start writes
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto* ppp = isA<processorLduInterface>(intf);
        const auto& faceCells = intf.faceCells();

        auto& sendBuf = sendBufs_[inti];

        sendBuf.resize_nocopy(faceCells.size());
        forAll(faceCells, face)
        {
            sendBuf[face] = psiInternal[faceCells[face]];
        }

        requests.push_back(UPstream::Request());
        UOPstream::write
        (
            requests.back(),
            ppp->neighbProcNo(),
            sendBuf.cdata_bytes(),
            sendBuf.size_bytes(),
            ppp->tag()+70,          // random offset
            ppp->comm()
        );
    }
}


void Foam::distributedDILUPreconditioner::wait
(
    DynamicList<UPstream::Request>& requests,
    const bool cancel
) const
{
    if (cancel)
    {
        UPstream::cancelRequests(requests);
    }
    else
    {
        UPstream::waitRequests(requests);
    }
    requests.clear();
}


// Diagonal calculation
// ~~~~~~~~~~~~~~~~~~~~

void Foam::distributedDILUPreconditioner::addInterfaceDiag
(
    solveScalarField& rD,
    const label inti,
    const Field<solveScalar>& recvBuf
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();
    const auto& interfaceIntCoeffs = solver_.interfaceIntCoeffs();

    const auto& intf = interfaces[inti].interface();
    // TBD: do not use patch faceCells but passed-in addressing?
    const auto& faceCells = intf.faceCells();
    const auto& bc = interfaceBouCoeffs[inti];
    const auto& ic = interfaceIntCoeffs[inti];

    forAll(recvBuf, face)
    {
        // Note:interfaceBouCoeffs, interfaceIntCoeffs on the receiving
        //      side are negated
        rD[faceCells[face]] -= bc[face]*ic[face]/recvBuf[face];
    }
}


void Foam::distributedDILUPreconditioner::forwardInternalDiag
(
    solveScalarField& rD,
    const label colouri
) const
{
    // Add forward constributions to diagonal

    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix.lower().begin();

    const label nFaces = matrix.upper().size();
    if (cellColourPtr_)
    {
        const auto& cellColour = *cellColourPtr_;
        for (label face=0; face<nFaces; face++)
        {
            const label cell = lPtr[face];
            if (cellColour[cell] == colouri)
            {
                rD[uPtr[face]] -= upperPtr[face]*lowerPtr[face]/rD[cell];
            }
        }
    }
    else
    {
        for (label face=0; face<nFaces; face++)
        {
            rD[uPtr[face]] -= upperPtr[face]*lowerPtr[face]/rD[lPtr[face]];
        }
    }
}


void Foam::distributedDILUPreconditioner::addInterface
(
    solveScalarField& wA,
    const label inti,
    const Field<solveScalar>& recvBuf
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();

    const auto& intf = interfaces[inti].interface();
    const auto& faceCells = intf.faceCells();
    const auto& bc = interfaceBouCoeffs[inti];
    const solveScalar* __restrict__ rDPtr = rD_.begin();
    solveScalar* __restrict__ wAPtr = wA.begin();

    forAll(recvBuf, face)
    {
        // Note: interfaceBouCoeffs is -upperPtr
        const label cell = faceCells[face];
        wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
    }
}


void Foam::distributedDILUPreconditioner::forwardInternal
(
    solveScalarField& wA,
    const label colouri
) const
{
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const label* const __restrict__ losortPtr = lduAddr.losortAddr().begin();
    const scalar* const __restrict__ lowerPtr = matrix.lower().begin();

    const label nFaces = matrix.upper().size();
    if (cellColourPtr_)
    {
        const auto& cellColour = *cellColourPtr_;
        for (label face=0; face<nFaces; face++)
        {
            const label sface = losortPtr[face];
            const label cell = lPtr[sface];
            if (cellColour[cell] == colouri)
            {
                wAPtr[uPtr[sface]] -=
                    rDPtr[uPtr[sface]]*lowerPtr[sface]*wAPtr[cell];
            }
        }
    }
    else
    {
        for (label face=0; face<nFaces; face++)
        {
            const label sface = losortPtr[face];
            wAPtr[uPtr[sface]] -=
                rDPtr[uPtr[sface]]*lowerPtr[sface]*wAPtr[lPtr[sface]];
        }
    }
}


void Foam::distributedDILUPreconditioner::backwardInternal
(
    solveScalarField& wA,
    const label colouri
) const
{
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    const label nFacesM1 = matrix.upper().size() - 1;

    if (cellColourPtr_)
    {
        const auto& cellColour = *cellColourPtr_;
        for (label face=nFacesM1; face>=0; face--)
        {
            const label cell = uPtr[face];
            if (cellColour[cell] == colouri)
            {
                // Note: lower cell guaranteed in same colour
                wAPtr[lPtr[face]] -=
                    rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[cell];
            }
        }
    }
    else
    {
        for (label face=nFacesM1; face>=0; face--)
        {
            wAPtr[lPtr[face]] -=
                rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
        }
    }
}


void Foam::distributedDILUPreconditioner::calcReciprocalD
(
    solveScalarField& rD
) const
{
    const auto& matrix = solver_.matrix();

    // Make sure no outstanding receives
    wait(lowerRecvRequests_);

    // Start reads (into recvBufs_)
    receive(lowerNbrs_, lowerRecvRequests_);

    // Start with diagonal
    const scalarField& diag = matrix.diag();
    rD.resize_nocopy(diag.size());
    std::copy(diag.begin(), diag.end(), rD.begin());


    // Subtract coupled contributions
    {
        // Wait for finish. Received result in recvBufs
        wait(lowerRecvRequests_);

        for (const label inti : lowerNbrs_)
        {
            addInterfaceDiag(rD, inti, recvBufs_[inti]);
        }
    }


    // - divide the internal mesh into domains/colour, similar to domain
    //   decomposition
    // - we could use subsetted bits of the mesh or just loop over only
    //   the cells of the current domain
    // - a domain only uses boundary values of lower numbered domains
    // - this also limits the interfaces that need to be evaluated since
    //   we assume an interface only changes local faceCells and these
    //   are all of the same colour

    for (label colouri = 0; colouri < nColours_; colouri++)
    {
        if (cellColourPtr_) // && colourBufs_.set(colouri))
        {
            for (const label inti : lowerGlobalRecv_[colouri])
            {
                addInterfaceDiag(rD, inti, colourBufs_[colouri][inti]);
            }
        }

        forwardInternalDiag(rD, colouri);

        // Store effect of exchanging rD to higher interfaces in colourBufs_
        if (cellColourPtr_)
        {
            sendGlobal(higherGlobalSend_[colouri], rD, higherColour_[colouri]);
        }
    }

    // Make sure no outstanding sends
    wait(higherSendRequests_);

    // Start writes of rD (using sendBufs)
    send(higherNbrs_, rD, higherSendRequests_);

    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rD[cell] = 1.0/rD[cell];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedDILUPreconditioner::distributedDILUPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& dict
)
:
    lduMatrix::preconditioner(sol),
    coupled_(dict.getOrDefault<bool>("coupled", true, keyType::LITERAL)),
    rD_(sol.matrix().diag().size())
{
    const lduMesh& mesh = sol.matrix().mesh();
    const auto& interfaces = sol.interfaces();
    const auto& interfaceBouCoeffs = sol.interfaceBouCoeffs();

    // Allocate buffers
    // ~~~~~~~~~~~~~~~~

    sendBufs_.setSize(interfaces.size());
    recvBufs_.setSize(interfaces.size());
    forAll(interfaceBouCoeffs, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& bc = interfaceBouCoeffs[inti];
            sendBufs_.set(inti, new solveScalarField(bc.size(), Zero));
            recvBufs_.set(inti, new solveScalarField(bc.size(), Zero));
        }
    }


    // Determine processor colouring
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //const processorColour& colours = processorColour::New(mesh);
    //const label myColour = colours[Pstream::myProcNo()];
    if (&mesh != meshPtr_)
    {
        procColoursPtr_.reset(new labelList(Pstream::nProcs(mesh.comm())));
        const label nProcColours =
            processorColour::colour(mesh, procColoursPtr_());

        if (debug)
        {
            Info<< typeName << " : calculated processor colours:"
                << nProcColours << endl;
        }

        meshPtr_ = &mesh;
    }
    const labelList& procColours = procColoursPtr_();

    const label myColour = procColours[Pstream::myProcNo(mesh.comm())];

    bool haveGlobalCoupled = false;
    bool haveDistributedAMI = false;
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& intf = interfaces[inti].interface();
            const auto* ppp = isA<processorLduInterface>(intf);
            if (ppp)
            {
                const label nbrColour = procColours[ppp->neighbProcNo()];
                //const label nbrColour = ppp->neighbProcNo();
                if (nbrColour < myColour)
                {
                    lowerNbrs_.append(inti);
                }
                else if (nbrColour > myColour)
                {
                    higherNbrs_.append(inti);
                }
                else
                {
                    WarningInFunction
                        << "weird processorLduInterface"
                        << " from " << ppp->myProcNo()
                        << " to " << ppp->neighbProcNo()
                        << endl;
                }
            }
            else //if (isA<cyclicAMILduInterface>(intf))
            {
                haveGlobalCoupled = true;

                const auto* AMIpp = isA<cyclicAMILduInterface>(intf);
                if (AMIpp)
                {
                    //const auto& AMI =
                    //(
                    //    AMIpp->owner()
                    //  ? AMIpp->AMI()
                    //  : AMIpp->neighbPatch().AMI()
                    //);
                    //haveDistributedAMI = AMI.distributed();

                    if (debug)
                    {
                        Info<< typeName << " : interface " << inti
                            << " of type " << intf.type()
                            << " is distributed:" << haveDistributedAMI << endl;
                    }
                }
            }
        }
    }



    if (debug)
    {
        Info<< typeName
            << " : lower proc interfaces:" << flatOutput(lowerNbrs_)
            << " higher proc interfaces:" << flatOutput(higherNbrs_)
            << endl;
    }


    // Determine local colouring/zoneing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off with all cells colour 0
    nColours_ = 1;
    cellColourPtr_.reset();
    if (coupled_ && haveGlobalCoupled)
    {
        labelList patchToColour;
        cellColourPtr_.reset(new labelList(mesh.lduAddr().size()));
        //if (haveDistributedAMI)
        //{
        //    nColours_ = 1;
        //    *cellColourPtr_ = 0;
        //    patchToColour = labelList(interfaces.size(), 0);
        //}
        //else
        {
            nColours_ = processorColour::cellColour
            (
                mesh,
                cellColourPtr_(),
                patchToColour
            );
        }

        lowerGlobalRecv_.resize_nocopy(nColours_);
        lowerGlobalSend_.resize_nocopy(nColours_);
        lowerColour_.setSize(nColours_, -1);
        higherGlobalRecv_.resize_nocopy(nColours_);
        higherGlobalSend_.resize_nocopy(nColours_);
        higherColour_.setSize(nColours_, -1);
        colourBufs_.setSize(nColours_);

        forAll(interfaces, inti)
        {
            if (interfaces.set(inti))
            {
                const auto& intf = interfaces[inti].interface();

                label nbrInti = -1;
                const auto* AMIpp = isA<cyclicAMILduInterface>(intf);
                if (AMIpp)
                {
                    nbrInti = AMIpp->neighbPatchID();
                }
                const auto* cycpp = isA<cyclicLduInterface>(intf);
                if (cycpp)
                {
                    nbrInti = cycpp->neighbPatchID();
                }

                if (nbrInti != -1)
                {
                    const label colouri = patchToColour[inti];
                    const label nbrColouri = patchToColour[nbrInti];

                    if
                    (
                        (haveDistributedAMI && (inti < nbrInti))
                     || (!haveDistributedAMI && (colouri < nbrColouri))
                    )
                    {
                        // Send to higher
                        higherGlobalSend_[colouri].append(nbrInti);
                        higherColour_[colouri] = nbrColouri;
                        // Receive from higher
                        higherGlobalRecv_[colouri].append(inti);
                    }
                    else
                    {
                        // Send to lower
                        lowerGlobalSend_[colouri].append(nbrInti);
                        lowerColour_[colouri] = nbrColouri;
                        // Receive from lower
                        lowerGlobalRecv_[colouri].append(inti);
                    }
                }
            }
        }
    }


    if (debug)
    {
        Info<< typeName << " : local colours:" << nColours_ << endl;
        Info<< incrIndent;
        forAll(lowerGlobalRecv_, colouri)
        {
            const auto& lowerRecv = lowerGlobalRecv_[colouri];
            if (lowerRecv.size())
            {
                const auto& lowerSend = lowerGlobalSend_[colouri];
                const auto& lowerCol = lowerColour_[colouri];
                Info<< indent
                    << "Lower global interfaces for colour:" << colouri
                    << " receiving from:" << flatOutput(lowerRecv)
                    << " sending to:" << flatOutput(lowerSend)
                    << " with colour:" << lowerCol << endl;
            }
        }
        forAll(higherGlobalRecv_, colouri)
        {
            const auto& higherRecv = higherGlobalRecv_[colouri];
            if (higherRecv.size())
            {
                const auto& higherSend = higherGlobalSend_[colouri];
                const auto& higherCol = higherColour_[colouri];
                Info<< indent
                    << "Higher global interfaces for colour:" << colouri
                    << " receiving from:" << flatOutput(higherRecv)
                    << " sending to:" << flatOutput(higherSend)
                    << " with colour:" << higherCol << endl;
            }
        }
    }

    calcReciprocalD(rD_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedDILUPreconditioner::~distributedDILUPreconditioner()
{
    DebugPout<< "~distributedDILUPreconditioner()" << endl;

    // Wait on all requests before storage is being taken down

    wait(lowerSendRequests_);
    wait(lowerRecvRequests_);
    wait(higherSendRequests_);
    wait(higherRecvRequests_);

    // TBD: cancel/ignore outstanding messages
    //wait(lowerSendRequests_, true);
    //wait(lowerRecvRequests_, true);
    //wait(higherSendRequests_, true);
    //wait(higherRecvRequests_, true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributedDILUPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rAPtr = rA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label nCells = wA.size();


    // Forward sweep
    // ~~~~~~~~~~~~~

    // Make sure no receives are still in flight (should not happen)
    wait(lowerRecvRequests_);

    // Start reads (into recvBufs)
    receive(lowerNbrs_, lowerRecvRequests_);

    // Initialise 'internal' cells
    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    // Do 'halo' contributions from lower numbered procs
    {
        // Wait for finish. Received result in recvBufs
        wait(lowerRecvRequests_);

        for (const label inti : lowerNbrs_)
        {
            addInterface(wA, inti, recvBufs_[inti]);
        }

        for (label colouri = 0; colouri < nColours_; colouri++)
        {
            // Do non-processor boundaries for this colour
            if (cellColourPtr_)  // && colourBufs_.set(colouri))
            {
                for (const label inti : lowerGlobalRecv_[colouri])
                {
                    addInterface(wA, inti, colourBufs_[colouri][inti]);
                }
            }

            forwardInternal(wA, colouri);

            // Store effect of exchanging rD to higher interfaces
            // in colourBufs_
            if (cellColourPtr_)
            {
                sendGlobal
                (
                    higherGlobalSend_[colouri],
                    wA,
                    higherColour_[colouri]
                );
            }
        }
    }

    // Make sure no outstanding sends from previous iteration
    wait(higherSendRequests_);

    // Start writes of wA (using sendBufs)
    send(higherNbrs_, wA, higherSendRequests_);


    // Backward sweep
    // ~~~~~~~~~~~~~~

    // Make sure no outstanding receives
    wait(higherRecvRequests_);

    // Start receives
    receive(higherNbrs_, higherRecvRequests_);

    {
        // Wait until receives finished
        wait(higherRecvRequests_);

        for (const label inti : higherNbrs_)
        {
            addInterface(wA, inti, recvBufs_[inti]);
        }

        for (label colouri = nColours_-1; colouri >= 0; colouri--)
        {
            // Do non-processor boundaries for this colour
            if (cellColourPtr_)  // && colourBufs_.set(colouri))
            {
                for (const label inti : higherGlobalRecv_[colouri])
                {
                    addInterface(wA, inti, colourBufs_[colouri][inti]);
                }
            }

            backwardInternal(wA, colouri);

            // Store effect of exchanging rD to higher interfaces
            // in colourBufs_
            if (cellColourPtr_)
            {
                sendGlobal
                (
                    lowerGlobalSend_[colouri],
                    wA,
                    lowerColour_[colouri]
                );
            }
        }
    }

    // Make sure no outstanding sends
    wait(lowerSendRequests_);

    // Start writes of wA (using sendBufs)
    send(lowerNbrs_, wA, lowerSendRequests_);
}


void Foam::distributedDILUPreconditioner::setFinished
(
    const solverPerformance& s
) const
{
    DebugPout<< "setFinished fieldName:" << s.fieldName() << endl;

    // Wait on all requests before storage is being taken down
    // (could rely on construction order?)
    wait(lowerSendRequests_);
    wait(lowerRecvRequests_);
    wait(higherSendRequests_);
    wait(higherRecvRequests_);

    // TBD: cancel/ignore outstanding messages
    //wait(lowerSendRequests_, true);
    //wait(lowerRecvRequests_, true);
    //wait(higherSendRequests_, true);
    //wait(higherRecvRequests_, true);
}


// ************************************************************************* //
