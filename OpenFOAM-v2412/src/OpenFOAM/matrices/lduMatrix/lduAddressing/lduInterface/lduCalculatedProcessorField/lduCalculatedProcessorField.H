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

Class
    Foam::lduCalculatedProcessorField

Group
    grpGenericBoundaryConditions

Description
    A lduProcessorField type bypassing coupledFvPatchField

    Used to add updateInterfaceMatrix capabilities to a lduMatrix
    which is fully uncoupled from the fvMesh.

    Its functionality is purely to init and update the processor interfaces.

SourceFiles
    lduCalculatedProcessorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduCalculatedProcessorField_H
#define Foam_lduCalculatedProcessorField_H

#include "lduPrimitiveProcessorInterface.H"
#include "processorLduInterfaceField.H"
#include "LduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class lduCalculatedProcessorField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lduCalculatedProcessorField
:
    public LduInterfaceField<Type>,
    public processorLduInterfaceField
{
protected:

    // Protected Data

        //- Local reference cast into the interface
         const lduPrimitiveProcessorInterface& procInterface_;


        // Sending and receiving

            //- Current (non-blocking) send request
            mutable label sendRequest_;

            //- Current (non-blocking) recv request
            mutable label recvRequest_;

            //- Send buffer
            mutable Field<Type> sendBuf_;

            //- Receive buffer
            mutable Field<Type> recvBuf_;

            //- Scalar send buffer
            mutable solveScalarField scalarSendBuf_;

            //- Scalar recv buffer
            mutable solveScalarField scalarRecvBuf_;


    // Protected Member Functions

        void addToInternalField
        (
            solveScalarField& result,
            const bool add,
            const scalarField& coeffs,
            const solveScalarField& vals
        ) const;

        //- Receive and send requests have both completed
        virtual bool all_ready() const;


public:

    //- Runtime type information
    ClassName("lduCalculatedProcessorField");


    // Constructors

        //- Construct from ldu interface
        explicit lduCalculatedProcessorField
        (
            const lduInterface& interface
        );

        //- Copy construct
        lduCalculatedProcessorField
        (
            const lduCalculatedProcessorField<Type>&
        );


    //- Destructor
    virtual ~lduCalculatedProcessorField() = default;


    // Member Functions

    // Evaluation

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Initialise neighbour matrix update
        virtual void initInterfaceMatrixUpdate
        (
            solveScalarField& result,
            const bool add,
            const lduAddressing& lduAddr,
            const label patchId,
            const solveScalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            solveScalarField& result,
            const bool add,
            const lduAddressing& lduAddr,
            const label patchId,
            const solveScalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        //- Initialise neighbour matrix update
        virtual void initInterfaceMatrixUpdate
        (
            Field<scalar>& result,
            const bool add,
            const lduAddressing& lduAddr,
            const label patchId,
            const Field<scalar>& psiInternal,
            const scalarField& coeffs,
            const Pstream::commsTypes commsType
        ) const
        {
            NotImplemented;
        }

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            Field<scalar>& result,
            const bool add,
            const lduAddressing& lduAddr,
            const label patchId,
            const Field<scalar>& psiInternal,
            const scalarField& coeffs,
            const Pstream::commsTypes commsType
        ) const
        {
            NotImplemented;
        }


    // Processor coupled interface functions

        //- Return communicator used for communication
        virtual label comm() const
        {
            return procInterface_.comm();
        }

        //- Return processor number
        virtual int myProcNo() const
        {
            return procInterface_.myProcNo();
        }

        //- Return neighbour processor number
        virtual int neighbProcNo() const
        {
            return procInterface_.myProcNo();
        }

        //- Is the transform required
        virtual bool doTransform() const
        {
            return false;
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return procInterface_.forwardT();
        }

        //- Return rank of component for transform
        virtual int rank() const
        {
            return pTraits<Type>::rank;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduCalculatedProcessorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
