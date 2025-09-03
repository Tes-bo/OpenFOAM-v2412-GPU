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
    Foam::cyclicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries, whereby communication between the patches is performed using
    an arbitrary mesh interface (AMI) interpolation.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicAMI;
        value           <initial value>;
        neighbourValue  <initial value of neighbour patch cells>;
    }
    \endverbatim

Note
    The outer boundary of the patch pairs must be similar, i.e. if the owner
    patch is transformed to the neighbour patch, the outer perimeter of each
    patch should be identical (or very similar).

    The \c neighbourValue is only needed when running distributed,
    i.e. the neighbour cells are on a different processor from the owner cells.
    It guarantees consistent restart e.g. when doing a snGrad and avoids
    additional communication.

See also
    Foam::AMIInterpolation
    Foam::processorFvPatchField

SourceFiles
    cyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicAMIFvPatchField_H
#define Foam_cyclicAMIFvPatchField_H

#include "cyclicAMIFvPatch.H"
#include "coupledFvPatchField.H"
#include "cyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicAMIFvPatchField
:
    virtual public cyclicAMILduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const cyclicAMIFvPatch& cyclicAMIPatch_;


        // Sending and receiving (distributed AMI)

            //- Current range of send requests (non-blocking)
            mutable labelRange sendRequests_;

            //- Current range of recv requests (non-blocking)
            mutable labelRange recvRequests_;

            //- Send buffers
            mutable PtrList<List<Type>> sendBufs_;

            //- Receive buffers_
            mutable PtrList<List<Type>> recvBufs_;

            //- Scalar send buffers
            mutable PtrList<List<solveScalar>> scalarSendBufs_;

            //- Scalar receive buffers
            mutable PtrList<List<solveScalar>> scalarRecvBufs_;

            //- Neighbour coupled internal cell data
            mutable autoPtr<Field<Type>> patchNeighbourFieldPtr_;


    // Private Member Functions

        //- Return the AMI corresponding to the owner side
        const AMIPatchToPatchInterpolation& ownerAMI() const
        {
            return
            (
                cyclicAMIPatch_.owner()
              ? cyclicAMIPatch_.AMI()
              : cyclicAMIPatch_.neighbPatch().AMI()
            );
        }

        //- All receive/send requests have completed
        virtual bool all_ready() const;

        //- Use neighbour field caching
        static bool cacheNeighbourField();

        //- Return neighbour coupled internal cell data
        tmp<Field<Type>> patchNeighbourField(const Field<Type>&) const;

        //- Return neighbour side field given internal fields
        template<class Type2>
        tmp<Field<Type2>> neighbourSideField
        (
            const Field<Type2>&
        ) const;


        //- Return new matrix coeffs
        tmp<Field<scalar>> coeffs
        (
            fvMatrix<Type>& matrix,
            const Field<scalar>&,
            const label
        ) const;

        template<class Type2>
        void collectStencilData
        (
            const refPtr<mapDistribute>& mapPtr,
            const labelListList& stencil,
            const Type2& data,
            List<Type2>& expandedData
        );


public:

    //- Runtime type information
    TypeName(cyclicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicAMIFvPatchField onto a new patch
        cyclicAMIFvPatchField
        (
            const cyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicAMIFvPatchField(const cyclicAMIFvPatchField<Type>&);

        //- Construct as copy setting internal field reference
        cyclicAMIFvPatchField
        (
            const cyclicAMIFvPatchField<Type>&,
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


    // Member Functions

        //- Return local reference cast into the cyclic AMI patch
        const cyclicAMIFvPatch& cyclicAMIPatch() const noexcept
        {
            return cyclicAMIPatch_;
        }


    // Coupling

        //- Return true if coupled. Note that the underlying patch
        //  is not coupled() - the points don't align.
        virtual bool coupled() const { return cyclicAMIPatch_.coupled(); }

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Return neighbour coupled internal cell data
        virtual tmp<Field<Type>> patchNeighbourField() const;

        //- Return reference to neighbour patchField
        const cyclicAMIFvPatchField<Type>& neighbourPatchField() const;


    // Mapping Functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<Type>&,
            const labelList&
        );


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
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;


            //- Manipulate matrix
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& m,
                const label iMatrix,
                const direction cmpt
            );


        // Coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return (pTraits<Type>::rank && !cyclicAMIPatch_.parallel());
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicAMIPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicAMIPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;


        // Member Operators

            virtual void operator=(const fvPatchField<Type>&);
            virtual void operator==(const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
