/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "cyclicACMIGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterfaceField
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0),
    sendRequests_(),
    recvRequests_()
{
    const cyclicAMILduInterfaceField& p =
        refCast<const cyclicAMILduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank),
    sendRequests_(),
    recvRequests_()
{}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    Istream& is
)
:
    GAMGInterfaceField(GAMGCp, is),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(readBool(is)),
    rank_(readLabel(is)),
    sendRequests_(),
    recvRequests_()
{}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& local,
    const UPtrList<lduInterfaceField>& other
)
:
    GAMGInterfaceField(GAMGCp, local),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0),
    sendRequests_(),
    recvRequests_()
{
    const auto& p = refCast<const cyclicACMILduInterfaceField>(local);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicACMIGAMGInterfaceField::ready() const
{
    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
    )
    {
        recvRequests_.clear();

        if
        (
            UPstream::finishedRequests
            (
                sendRequests_.start(),
                sendRequests_.size()
            )
        )
        {
            sendRequests_.clear();
        }

        return true;
    }

    return false;
}


void Foam::cyclicACMIGAMGInterfaceField::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const auto& AMI =
    (
        cyclicACMIInterface_.owner()
      ? cyclicACMIInterface_.AMI()
      : cyclicACMIInterface_.neighbPatch().AMI()
    );

    if (AMI.distributed())
    {
        DebugPout<< "cyclicACMIFvPatchField::initInterfaceMatrixUpdate() :"
            << " interface:" << cyclicACMIInterface_.index()
            << " size:" << cyclicACMIInterface_.size()
            << " owner:" << cyclicACMIInterface_.owner()
            << " AMI distributed:" << AMI.distributed()
            << endl;

        // Start sending
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        // Get neighbouring field
        const labelList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIInterface_.neighbPatchID());

        solveScalarField pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        const auto& map =
        (
            cyclicACMIInterface_.owner()
          ? AMI.tgtMap()
          : AMI.srcMap()
        );

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIInterface_.index()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();

        // Insert send/receive requests (non-blocking). See e.g.
        // cyclicAMIPolyPatchTemplates.C
        const label oldWarnComm = UPstream::commWarn(AMI.comm());
        map.send
        (
            pnf,
            sendRequests_,
            scalarSendBufs_,
            recvRequests_,
            scalarRecvBufs_,
            19462+cyclicACMIInterface_.index()  // unique offset + patch index
        );
        UPstream::commWarn(oldWarnComm);
    }

    this->updatedMatrix(false);
}


void Foam::cyclicACMIGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    const auto& AMI =
    (
        cyclicACMIInterface_.owner()
      ? cyclicACMIInterface_.AMI()
      : cyclicACMIInterface_.neighbPatch().AMI()
    );

    DebugPout<< "cyclicACMIGAMGInterfaceField::updateInterfaceMatrix() :"
        << " interface:" << cyclicACMIInterface_.index()
        << " size:" << cyclicACMIInterface_.size()
        << " owner:" << cyclicACMIInterface_.owner()
        << " AMI distributed:" << AMI.distributed()
        << endl;


    if (AMI.distributed())
    {
        const auto& map =
        (
            cyclicACMIInterface_.owner()
          ? AMI.tgtMap()
          : AMI.srcMap()
        );

        // Receive (= copy) data from buffers into work. TBD: receive directly
        // into slices of work.
        solveScalarField work;
        map.receive
        (
            recvRequests_,
            scalarRecvBufs_,
            work,
            19462+cyclicACMIInterface_.index()  // unique offset + patch index
        );

        // Receive requests all handled by last function call
        recvRequests_.clear();

        solveScalarField pnf(faceCells.size(), Zero);
        AMI.weightedSum
        (
            cyclicACMIInterface_.owner(),
            work,
            pnf,               // result
            solveScalarField::null()
        );

        // Add result using coefficients
        this->addToInternalField(result, !add, faceCells, coeffs, pnf);
    }
    else
    {
        // Get neighbouring field
        const labelList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIInterface_.neighbPatchID());

        solveScalarField pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        if (cyclicACMIInterface_.owner())
        {
            pnf = AMI.interpolateToSource(pnf);
        }
        else
        {
            pnf = AMI.interpolateToTarget(pnf);
        }

        const labelUList& faceCells = lduAddr.patchAddr(patchId);

        this->addToInternalField(result, !add, faceCells, coeffs, pnf);
    }

    this->updatedMatrix(true);
}


void Foam::cyclicACMIGAMGInterfaceField::write(Ostream& os) const
{
    //GAMGInterfaceField::write(os);
    os  << token::SPACE << doTransform()
        << token::SPACE << rank();
}


// ************************************************************************* //
