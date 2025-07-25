/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::cyclicFaPatchField

Description
    Foam::cyclicFaPatchField

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicFaPatchField_H
#define Foam_cyclicFaPatchField_H

#include "coupledFaPatchField.H"
#include "cyclicLduInterfaceField.H"
#include "cyclicFaPatch.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cyclicFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFaPatchField
:
    virtual public cyclicLduInterfaceField,
    public coupledFaPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const cyclicFaPatch& cyclicPatch_;


    // Private Member Functions

        //- Return neighbour side field given internal fields
        template<class Type2>
        tmp<Field<Type2>> neighbourSideField
        (
            const Field<Type2>&
        ) const;


public:

    //- Runtime type information
    TypeName(cyclicFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary& dict,
            //! Evaluate (default: mandatory)
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Construct by mapping given cyclicFaPatchField onto a new patch
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&,
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


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicFaPatch& cyclicPatch() const
            {
                return cyclicPatch_;
            }


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type>> patchNeighbourField() const;


        // Coupled interface functionality

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


        //- Cyclic coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return (pTraits<Type>::rank && !cyclicPatch_.parallel());
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPatch_.reverseT();
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
    #include "cyclicFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
