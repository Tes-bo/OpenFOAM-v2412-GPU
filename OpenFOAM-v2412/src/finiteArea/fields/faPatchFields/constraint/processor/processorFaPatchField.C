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

#include "processorFaPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    coupledFaPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFaPatch>(p)),
    sendRequest_(-1),
    recvRequest_(-1)
{
    if (!isType<processorFaPatch>(this->patch()))
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
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    procPatch_(refCast<const processorFaPatch>(p, dict)),
    sendRequest_(-1),
    recvRequest_(-1)
{
    if (!isType<processorFaPatch>(p))
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
        this->extrapolateInternal();  // Zero-gradient patch values
    }
}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
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
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
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
bool Foam::processorFaPatchField<Type>::all_ready() const
{
    return UPstream::finishedRequestPair(recvRequest_, sendRequest_);
}


template<class Type>
bool Foam::processorFaPatchField<Type>::ready() const
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
Foam::processorFaPatchField<Type>::patchNeighbourField() const
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
void Foam::processorFaPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (UPstream::parRun())
    {
        this->patchInternalField(sendBuf_);

        if (commsType == UPstream::commsTypes::nonBlocking)
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
            procPatch_.send(commsType, sendBuf_);
        }
    }
}


template<class Type>
void Foam::processorFaPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (UPstream::parRun())
    {
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            // Fast path. Received into *this

            // Require receive data.
            // Only update the send request state.
            UPstream::waitRequest(recvRequest_); recvRequest_ = -1;
            if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
        }
        else
        {
            procPatch_.receive<Type>(commsType, *this);
        }

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorFaPatchField<Type>::snGrad() const
{
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, scalarSendBuf_);

    if (commsType == UPstream::commsTypes::nonBlocking)
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
        procPatch_.send(commsType, scalarSendBuf_);
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
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

    const labelUList& faceCells = this->patch().edgeFaces();

    if (commsType == UPstream::commsTypes::nonBlocking)
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
        procPatch_.receive(commsType, scalarRecvBuf_);
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
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, sendBuf_);

    if (commsType == UPstream::commsTypes::nonBlocking)
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
        procPatch_.send(commsType, sendBuf_);
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
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

    const labelUList& faceCells = this->patch().edgeFaces();

    if (commsType == UPstream::commsTypes::nonBlocking)
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
        procPatch_.receive(commsType, recvBuf_);
    }


    // Transform according to the transformation tensor
    transformCoupleField(recvBuf_);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, recvBuf_);

    this->updatedMatrix(true);
}


// ************************************************************************* //
