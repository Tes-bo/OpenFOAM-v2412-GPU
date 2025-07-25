/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "calculatedProcessorGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calculatedProcessorGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        calculatedProcessorGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        calculatedProcessorGAMGInterfaceField,
        lduInterfaceField
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        calculatedProcessorGAMGInterfaceField,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calculatedProcessorGAMGInterfaceField::
calculatedProcessorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    procInterface_(refCast<const calculatedProcessorGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0),
    sendRequest_(-1),
    recvRequest_(-1)
{
    const auto& p = refCast<const processorLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::calculatedProcessorGAMGInterfaceField::
calculatedProcessorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    procInterface_(refCast<const calculatedProcessorGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank),
    sendRequest_(-1),
    recvRequest_(-1)
{}


Foam::calculatedProcessorGAMGInterfaceField::
calculatedProcessorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    Istream& is
)
:
    GAMGInterfaceField(GAMGCp, is),
    procInterface_(refCast<const calculatedProcessorGAMGInterface>(GAMGCp)),
    doTransform_(readBool(is)),
    rank_(readLabel(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calculatedProcessorGAMGInterfaceField::initInterfaceMatrixUpdate
(
    solveScalarField&,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    procInterface_.interfaceInternalField(psiInternal, scalarSendBuf_);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !UPstream::floatTransfer
    )
    {
        // Fast path.
        scalarRecvBuf_.resize_nocopy(scalarSendBuf_.size());

        recvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            scalarRecvBuf_.data_bytes(),
            scalarRecvBuf_.size_bytes(),
            procInterface_.tag(),
            comm()
        );

        sendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            scalarSendBuf_.cdata_bytes(),
            scalarSendBuf_.size_bytes(),
            procInterface_.tag(),
            comm()
        );
    }
    else
    {
        procInterface_.compressedSend(commsType, scalarSendBuf_);
    }

    this->updatedMatrix(false);
}


void Foam::calculatedProcessorGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !UPstream::floatTransfer
    )
    {
        // Fast path: consume straight from receive buffer

        // Require receive data.
        // Only update the send request state.
        UPstream::waitRequest(recvRequest_); recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }
    else
    {
        scalarRecvBuf_.resize_nocopy(this->size());
        procInterface_.compressedReceive(commsType, scalarRecvBuf_);
    }


    // Transform according to the transformation tensor
    transformCoupleField(scalarRecvBuf_, cmpt);

    // Multiply the field by coefficients and add into the result
    addToInternalField(result, !add, faceCells, coeffs, scalarRecvBuf_);

    this->updatedMatrix(true);
}


void Foam::calculatedProcessorGAMGInterfaceField::write(Ostream& os) const
{
    //GAMGInterfaceField::write(os);
    os  << token::SPACE << doTransform()
        << token::SPACE << rank();
}


// ************************************************************************* //
