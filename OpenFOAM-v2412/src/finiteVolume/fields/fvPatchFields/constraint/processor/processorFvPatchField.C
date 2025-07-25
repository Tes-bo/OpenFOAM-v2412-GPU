/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    procPatch_(refCast<const processorFvPatch>(p, dict)),
    sendRequest_(-1),
    recvRequest_(-1)
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    // Use 'value' supplied, or set to internal field
    if (!this->readValueEntry(dict))
    {
        this->extrapolateInternal();
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << procPatch_.name()
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    sendRequest_(-1),
    recvRequest_(-1),
    sendBuf_(std::move(ptf.sendBuf_)),
    recvBuf_(std::move(ptf.recvBuf_)),
    scalarSendBuf_(std::move(ptf.scalarSendBuf_)),
    scalarRecvBuf_(std::move(ptf.scalarRecvBuf_))
{
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << procPatch_.name()
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    sendRequest_(-1),
    recvRequest_(-1)
{
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << procPatch_.name()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::processorFvPatchField<Type>::all_ready() const
{
    return UPstream::finishedRequestPair(recvRequest_, sendRequest_);
}


template<class Type>
bool Foam::processorFvPatchField<Type>::ready() const
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
Foam::tmp<Foam::Field<Type>>
Foam::processorFvPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "Outstanding request on patch " << procPatch_.name()
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (UPstream::parRun())
    {
        this->patchInternalField(sendBuf_);

        if
        (
            commsType == UPstream::commsTypes::nonBlocking
         && (std::is_integral<Type>::value || !UPstream::floatTransfer)
        )
        {
            if (!is_contiguous<Type>::value)
            {
                FatalErrorInFunction
                    << "Invalid for non-contiguous data types"
                    << abort(FatalError);
            }

            // Receive straight into *this
            this->resize_nocopy(sendBuf_.size());

            recvRequest_ = UPstream::nRequests();
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                this->data_bytes(),
                this->size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );

            sendRequest_ = UPstream::nRequests();
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                sendBuf_.cdata_bytes(),
                sendBuf_.size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        else
        {
            procPatch_.compressedSend(commsType, sendBuf_);
        }
    }
}


template<class Type>
void Foam::processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (UPstream::parRun())
    {
        if
        (
            commsType == UPstream::commsTypes::nonBlocking
         && (std::is_integral<Type>::value || !UPstream::floatTransfer)
        )
        {
            // Fast path: received into *this

            // Require receive data.
            // Only update the send request state.
            UPstream::waitRequest(recvRequest_); recvRequest_ = -1;
            if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
        }
        else
        {
            procPatch_.compressedReceive<Type>(commsType, *this);
        }

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::processorFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    return deltaCoeffs*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
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
    //this->patch().patchInternalField(psiInternal, scalarSendBuf_);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    scalarSendBuf_.resize_nocopy(this->patch().size());
    forAll(scalarSendBuf_, facei)
    {
        scalarSendBuf_[facei] = psiInternal[faceCells[facei]];
    }

    if
    (
        commsType == UPstream::commsTypes::nonBlocking
     && !UPstream::floatTransfer
    )
    {
        // Fast path.
        if (debug && !this->all_ready())
        {
            FatalErrorInFunction
                << "Outstanding request(s) on patch " << procPatch_.name()
                << abort(FatalError);
        }

        scalarRecvBuf_.resize_nocopy(scalarSendBuf_.size());

        recvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            scalarRecvBuf_.data_bytes(),
            scalarRecvBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        sendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            scalarSendBuf_.cdata_bytes(),
            scalarSendBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, scalarSendBuf_);
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
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
        commsType == UPstream::commsTypes::nonBlocking
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
        procPatch_.compressedReceive(commsType, scalarRecvBuf_);
    }


    if (pTraits<Type>::rank)
    {
        // Transform non-scalar data according to the transformation tensor
        transformCoupleField(scalarRecvBuf_, cmpt);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, scalarRecvBuf_);

    this->updatedMatrix(true);
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>&,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField&,
    const Pstream::commsTypes commsType
) const
{
    sendBuf_.resize_nocopy(this->patch().size());

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    forAll(sendBuf_, facei)
    {
        sendBuf_[facei] = psiInternal[faceCells[facei]];
    }

    if
    (
        commsType == UPstream::commsTypes::nonBlocking
     && (std::is_integral<Type>::value || !UPstream::floatTransfer)
    )
    {
        // Fast path.
        if (debug && !this->all_ready())
        {
            FatalErrorInFunction
                << "Outstanding request(s) on patch " << procPatch_.name()
                << abort(FatalError);
        }

        recvBuf_.resize_nocopy(sendBuf_.size());

        recvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            recvBuf_.data_bytes(),
            recvBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        sendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, sendBuf_);
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>&,
    const scalarField& coeffs,
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
        commsType == UPstream::commsTypes::nonBlocking
     && (std::is_integral<Type>::value || !UPstream::floatTransfer)
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
        recvBuf_.resize_nocopy(this->size());
        procPatch_.compressedReceive(commsType, recvBuf_);
    }


    // Transform according to the transformation tensor
    transformCoupleField(recvBuf_);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, recvBuf_);

    this->updatedMatrix(true);
}


// ************************************************************************* //
