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

Class
    Foam::calculatedProcessorFvPatchField

Group
    grpGenericBoundaryConditions

Description
    A processorFvPatchField type bypassing fvPatch

    Used to temporarily add updateInterfaceMatrix capabilities to a matrix
    during overset solving. Supplies:
    - patchNeighbourField functionality (cached in *this as per
      processorFvPatchField)
    - initEvaluate/evaluate: caching of patchNeighbourField (see above)
    - initInterfaceMatrixUpdate etc: adding of neighbouring data

SourceFiles
    calculatedProcessorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_calculatedProcessorFvPatchField_H
#define Foam_calculatedProcessorFvPatchField_H

#include "lduPrimitiveProcessorInterface.H"
#include "coupledFvPatchField.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class calculatedProcessorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class calculatedProcessorFvPatchField
:
    public processorLduInterfaceField,
    public coupledFvPatchField<Type>
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
    TypeName("calculatedProcessor");


    // Constructors

        //- Construct from patch and internal field
        calculatedProcessorFvPatchField
        (
            const lduInterface& interface,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct as copy
        calculatedProcessorFvPatchField
        (
            const calculatedProcessorFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        calculatedProcessorFvPatchField
        (
            const calculatedProcessorFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~calculatedProcessorFvPatchField() = default;


    // Member Functions

    // Coupling

        //- The patch field is coupled if running in parallel
        virtual bool coupled() const { return UPstream::parRun(); }

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Return neighbour field of internal field
        virtual tmp<Field<Type>> patchNeighbourField() const;


    // Evaluation

        //- Initialise the evaluation of the patch field
        virtual void initEvaluate(const Pstream::commsTypes commsType);

        //- Evaluate the patch field
        virtual void evaluate(const Pstream::commsTypes commsType);


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
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
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
    #include "calculatedProcessorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
