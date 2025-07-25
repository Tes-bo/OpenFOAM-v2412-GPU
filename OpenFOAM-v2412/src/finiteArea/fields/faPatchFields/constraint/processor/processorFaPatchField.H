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

Class
    Foam::processorFaPatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    processorFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorFaPatchField_H
#define Foam_processorFaPatchField_H

#include "coupledFaPatchField.H"
#include "processorLduInterfaceField.H"
#include "processorFaPatch.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFaPatchField
:
    public processorLduInterfaceField,
    public coupledFaPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const processorFaPatch& procPatch_;


        // Sending and receiving

            //- Current (non-blocking) send request
            mutable label sendRequest_;

            //- Current (non-blocking) recv request
            mutable label recvRequest_;

            //- Send buffer.
            mutable Field<Type> sendBuf_;

            //- Receive buffer.
            mutable Field<Type> recvBuf_;

            //- Scalar send buffer
            mutable solveScalarField scalarSendBuf_;

            //- Scalar recv buffer
            mutable solveScalarField scalarRecvBuf_;


    // Private Member Functions

        //- Receive and send requests have both completed
        virtual bool all_ready() const;


public:

    //- Runtime type information
    TypeName(processorFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorFaPatchField onto a new patch
        processorFaPatchField
        (
            const processorFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        processorFaPatchField(const processorFaPatchField<Type>&);

        //- Construct as copy setting internal field reference
        processorFaPatchField
        (
            const processorFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return faPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    ~processorFaPatchField() = default;


    // Member Functions

    // Coupling

        //- The patch field is coupled if running in parallel
        virtual bool coupled() const { return UPstream::parRun(); }

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Return neighbour field given internal field
        tmp<Field<Type>> patchNeighbourField() const;


    // Evaluation

        //- Initialise the evaluation of the patch field
        virtual void initEvaluate(const Pstream::commsTypes commsType);

        //- Evaluate the patch field
        virtual void evaluate(const Pstream::commsTypes commsType);

        //- Initialise the evaluation of the patch field after a local
        //  operation. Dummy since operating on a copy
        virtual void initEvaluateLocal
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::buffered
        )
        {}

        //- Evaluate the patch field after a local operation (e.g. *=).
        //  Dummy since operating on a copy
        virtual void evaluateLocal
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::buffered
        )
        {}

        //- Return patch-normal gradient
        virtual tmp<Field<Type>> snGrad() const;


        // Coupled interface functionality

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
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>&,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>&,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;


    // Processor coupled interface functions

        //- Return communicator used for communication
        virtual label comm() const
        {
            return procPatch_.comm();
        }

        //- Return processor number
        virtual int myProcNo() const
        {
            return procPatch_.myProcNo();
        }

        //- Return neighbour processor number
        virtual int neighbProcNo() const
        {
            return procPatch_.neighbProcNo();
        }

        //- Does the patch field perform the transformation
        virtual bool doTransform() const
        {
            return (pTraits<Type>::rank && !procPatch_.parallel());
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return procPatch_.forwardT();
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
    #include "processorFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
