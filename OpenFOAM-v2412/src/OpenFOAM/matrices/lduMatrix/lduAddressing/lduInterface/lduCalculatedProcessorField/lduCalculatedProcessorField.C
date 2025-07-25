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

#include "lduCalculatedProcessorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduInterface& interface
)
:
    LduInterfaceField<Type>(interface),
    procInterface_(refCast<const lduPrimitiveProcessorInterface>(interface)),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduCalculatedProcessorField<Type>& ptf
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf)),
    procInterface_(ptf.procInterface_),
    sendRequest_(-1),
    recvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::lduCalculatedProcessorField<Type>::all_ready() const
{
    return UPstream::finishedRequestPair(recvRequest_, sendRequest_);
}


template<class Type>
bool Foam::lduCalculatedProcessorField<Type>::ready() const
{
    const bool ok = UPstream::finishedRequest(recvRequest_);
    if (ok)
    {
        recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }
    return ok;
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initInterfaceMatrixUpdate
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
    if (!this->all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on interface "
            //<< procInterface_.name()
            << abort(FatalError);
    }

    // Bypass patchInternalField since uses fvPatch addressing
    const labelList& fc = lduAddr.patchAddr(patchId);

    scalarSendBuf_.resize_nocopy(fc.size());
    forAll(fc, i)
    {
        scalarSendBuf_[i] = psiInternal[fc[i]];
    }

    scalarRecvBuf_.resize_nocopy(scalarSendBuf_.size());

    recvRequest_ = UPstream::nRequests();
    UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarRecvBuf_.data_bytes(),
        scalarRecvBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    sendRequest_ = UPstream::nRequests();
    UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarSendBuf_.cdata_bytes(),
        scalarSendBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    this->updatedMatrix(false);
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::addToInternalField
(
    solveScalarField& result,
    const bool add,
    const scalarField& coeffs,
    const solveScalarField& vals
) const
{
    const labelUList& faceCells = this->procInterface_.faceCells();

    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*vals[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*vals[elemI];
        }
    }
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::updateInterfaceMatrix
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
    if (this->updatedMatrix())
    {
        return;
    }

    {
        // Require receive data.
        // Only update the send request state.
        UPstream::waitRequest(recvRequest_); recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }

    // Consume straight from receive buffer. Note use of our own
    // helper to avoid using fvPatch addressing
    addToInternalField(result, !add, coeffs, scalarRecvBuf_);

    this->updatedMatrix(true);
}


// ************************************************************************* //
