/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::processorLduInterface

Description
    An abstract base class for processor coupled interfaces.

SourceFiles
    processorLduInterface.C
    processorLduInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorLduInterface_H
#define Foam_processorLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorLduInterface Declaration
\*---------------------------------------------------------------------------*/

class processorLduInterface
{
    // Private Data

        //- Send buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> byteSendBuf_;

        //- Receive buffer.
        //  Only sized and used when compressed or non-blocking comms used.
        mutable List<char> byteRecvBuf_;


    // Private Member Functions

        //- Increase buffer size if required
        //  Uses the nocopy variant since it will be overwritten
        static void resizeBuf(List<char>& buf, const label len)
        {
            if (buf.size() < len)
            {
                buf.resize_nocopy(len);
            }
        }


public:

    //- Runtime type information
    TypeNameNoDebug("processorLduInterface");


    // Constructors

        //- Default construct
        processorLduInterface() noexcept = default;


    //- Destructor
    virtual ~processorLduInterface() = default;


    // Member Functions

        // Access

            //- Return communicator used for parallel communication
            virtual label comm() const = 0;

            //- Return processor number (rank in communicator)
            virtual int myProcNo() const = 0;

            //- Return neighbour processor number (rank in communicator)
            virtual int neighbProcNo() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return message tag used for sending
            virtual int tag() const = 0;


        // Transfer Functions

            //- Raw send function
            template<class Type>
            void send
            (
                const UPstream::commsTypes commsType,
                const UList<Type>& f
            ) const;

            //- Raw receive function
            template<class Type>
            void receive
            (
                const UPstream::commsTypes commsType,
                UList<Type>& f
            ) const;

            //- Raw receive function returning field
            template<class Type>
            tmp<Field<Type>> receive
            (
                const UPstream::commsTypes commsType,
                const label size
            ) const;


            //- Raw send function with data compression
            template<class Type>
            void compressedSend
            (
                const UPstream::commsTypes commsType,
                const UList<Type>& f
            ) const;

            //- Raw receive function with data compression
            template<class Type>
            void compressedReceive
            (
                const UPstream::commsTypes commsType,
                UList<Type>& f
            ) const;

            //- Raw receive function with data compression returning field
            template<class Type>
            tmp<Field<Type>> compressedReceive
            (
                const UPstream::commsTypes commsType,
                const label size
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorLduInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
