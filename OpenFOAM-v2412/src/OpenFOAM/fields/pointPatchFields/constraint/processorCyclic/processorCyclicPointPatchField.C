/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "processorCyclicPointPatchField.H"
#include "transformField.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    procPatch_(refCast<const processorCyclicPointPatch>(p))
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorCyclicPointPatch>(p, dict))
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const processorCyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorCyclicPointPatch>(ptf.patch()))
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const processorCyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorCyclicPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorCyclicPointPatchField<Type>::initSwapAddSeparated
(
    const Pstream::commsTypes commsType,
    Field<Type>& pField
) const
{
    if (UPstream::parRun())
    {
        // Get internal field into correct order for opposite side. Note use
        // of member data buffer since using non-blocking. Could be optimised
        // out if not using non-blocking...
        this->patchInternalField
        (
            pField,
            procPatch_.reverseMeshPoints(),
            sendBuf_
        );

        if (commsType == Pstream::commsTypes::nonBlocking)
        {
            recvBuf_.resize_nocopy(sendBuf_.size());

            UIPstream::read
            (
                commsType,
                procPatch_.neighbProcNo(),
                recvBuf_.data_bytes(),
                recvBuf_.size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        UOPstream::write
        (
            commsType,
            procPatch_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
}


template<class Type>
void Foam::processorCyclicPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes commsType,
    Field<Type>& pField
) const
{
    if (UPstream::parRun())
    {
        // If nonblocking, data is already in receive buffer

        if (commsType != Pstream::commsTypes::nonBlocking)
        {
            recvBuf_.resize_nocopy(this->size());
            UIPstream::read
            (
                commsType,
                procPatch_.neighbProcNo(),
                recvBuf_.data_bytes(),
                recvBuf_.size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }

        if (doTransform())
        {
            const processorCyclicPolyPatch& ppp =
                procPatch_.procCyclicPolyPatch();
            const tensor& forwardT = ppp.forwardT()[0];

            transform(recvBuf_, forwardT, recvBuf_);
        }

        // All points are separated
        this->addToInternalField(pField, recvBuf_);
    }
}


// ************************************************************************* //
