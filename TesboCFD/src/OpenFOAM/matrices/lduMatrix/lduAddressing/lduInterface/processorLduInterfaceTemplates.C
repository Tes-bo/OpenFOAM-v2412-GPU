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

#include "processorLduInterface.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions * * *  * * * * * * * * * * //

template<class Type>
void Foam::processorLduInterface::send
(
    const UPstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    const label nBytes = f.byteSize();

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        UOPstream::write
        (
            commsType,
            neighbProcNo(),
            f.cdata_bytes(),
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        resizeBuf(byteSendBuf_, nBytes);
        std::memcpy
        (
            static_cast<void*>(byteSendBuf_.data()), f.cdata(), nBytes
        );

        resizeBuf(byteRecvBuf_, nBytes);

        if (!nBytes)
        {
            // Can skip empty messages
            return;
        }

        UIPstream::read
        (
            commsType,
            neighbProcNo(),
            byteRecvBuf_.data(),
            nBytes,
            tag(),
            comm()
        );

        UOPstream::write
        (
            commsType,
            neighbProcNo(),
            byteSendBuf_.cdata(),
            nBytes,
            tag(),
            comm()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
void Foam::processorLduInterface::receive
(
    const UPstream::commsTypes commsType,
    UList<Type>& f
) const
{
    const label nBytes = f.byteSize();

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        UIPstream::read
        (
            commsType,
            neighbProcNo(),
            f.data_bytes(),
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        std::memcpy
        (
            static_cast<void*>(f.data()), byteRecvBuf_.cdata(), nBytes
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::receive
(
    const UPstream::commsTypes commsType,
    const label size
) const
{
    auto tfld = tmp<Field<Type>>::New(size);
    receive(commsType, tfld.ref());
    return tfld;
}


template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const UPstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if
    (
        f.size()
     && UPstream::floatTransfer
     && (!std::is_integral<Type>::value && sizeof(scalar) != sizeof(float))
    )
    {
        static const label nCmpts =
        (
            // Placeholder value for unusable template instantiation
            std::is_integral<Type>::value
          ? 1
          : sizeof(Type)/sizeof(scalar)
        );
        const label nm1 = (f.size() - 1)*nCmpts;
        const label nBytes = f.size()*nCmpts*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.cdata());
        const scalar *slast = &sArray[nm1];
        resizeBuf(byteSendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(byteSendBuf_.data());

        for (label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        reinterpret_cast<Type&>(fArray[nm1]) = f.last();

        if
        (
            commsType == UPstream::commsTypes::buffered
         || commsType == UPstream::commsTypes::scheduled
        )
        {
            UOPstream::write
            (
                commsType,
                neighbProcNo(),
                byteSendBuf_.cdata(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == UPstream::commsTypes::nonBlocking)
        {
            resizeBuf(byteRecvBuf_, nBytes);

            UIPstream::read
            (
                commsType,
                neighbProcNo(),
                byteRecvBuf_.data(),
                nBytes,
                tag(),
                comm()
            );

            UOPstream::write
            (
                commsType,
                neighbProcNo(),
                byteSendBuf_.cdata(),
                nBytes,
                tag(),
                comm()
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }
    }
    else
    {
        this->send(commsType, f);
    }
}


template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const UPstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if
    (
        f.size()
     && UPstream::floatTransfer
     && (!std::is_integral<Type>::value && sizeof(scalar) != sizeof(float))
    )
    {
        static const label nCmpts =
        (
            // Placeholder value for unusable template instantiation
            std::is_integral<Type>::value
          ? 1
          : sizeof(Type)/sizeof(scalar)
        );
        const label nm1 = (f.size() - 1)*nCmpts;
        const label nBytes = f.size()*nCmpts*sizeof(float);

        if
        (
            commsType == UPstream::commsTypes::buffered
         || commsType == UPstream::commsTypes::scheduled
        )
        {
            resizeBuf(byteRecvBuf_, nBytes);

            UIPstream::read
            (
                commsType,
                neighbProcNo(),
                byteRecvBuf_.data(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }

        const float *fArray =
            reinterpret_cast<const float*>(byteRecvBuf_.cdata());
        f.last() = reinterpret_cast<const Type&>(fArray[nm1]);
        scalar *sArray = reinterpret_cast<scalar*>(f.data());
        const scalar *slast = &sArray[nm1];

        for (label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }
    }
    else
    {
        this->receive<Type>(commsType, f);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::compressedReceive
(
    const UPstream::commsTypes commsType,
    const label size
) const
{
    auto tfld = tmp<Field<Type>>::New(size);
    compressedReceive(commsType, tfld.ref());
    return tfld;
}


// ************************************************************************* //
