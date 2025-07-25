/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::mixedFaPatchField

Description
    This boundary condition provides a base class for 'mixed' type boundary
    conditions, i.e. conditions that mix fixed value and patch-normal gradient
    conditions.

Usage
    \table
        Property     | Description                          | Required | Default
        refValue     | fixed value                          | yes |
        refGradient  | patch normal gradient                | yes |
        valueFraction | value weighting (0-1)               | yes |
    \endtable

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixedFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mixedFaPatchField_H
#define Foam_mixedFaPatchField_H

#include "faPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixedFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedFaPatchField
:
    public faPatchField<Type>
{
    // Private Data

        //- Value field
        Field<Type> refValue_;

        //- Normal gradient field
        Field<Type> refGrad_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

protected:

        //- Read the "refValue", "refGradient" and "valueFraction" entries
        //- into their respective places.
        //  The reading can be optional (default), mandatory etc.
        //  If refValue is to be read, refGradient and valueFraction must
        //  also exist.
        //  \returns True on success
        bool readMixedEntries
        (
            const dictionary& dict,
            IOobjectOption::readOption readOpt = IOobjectOption::LAZY_READ
        );


public:

    //- Runtime type information
    TypeName("mixed");


    // Constructors

        //- Construct from patch and internal field
        mixedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field,
        //- initialise as zero-gradient
        mixedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Foam::zero
        );

        //- Construct from patch, internal field and dictionary
        mixedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&,
            //! The "refValue", "refGradient", "valueFraction" entries
            //! (default: mandatory)
            IOobjectOption::readOption requireMixed = IOobjectOption::MUST_READ
        );

        //- Construct by mapping the given mixedFaPatchField onto a new patch
        mixedFaPatchField
        (
            const mixedFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        mixedFaPatchField
        (
            const mixedFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        mixedFaPatchField
        (
            const mixedFaPatchField<Type>&,
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


    // Member Functions

        //- True: the patch field fixes a value.
        virtual bool fixesValue() const { return true; }


        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual Field<Type>& refGrad()
            {
                return refGrad_;
            }

            virtual const Field<Type>& refGrad() const
            {
                return refGrad_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const faPatchField<Type>&) {}
        virtual void operator+=(const faPatchField<Type>&) {}
        virtual void operator-=(const faPatchField<Type>&) {}
        virtual void operator*=(const faPatchField<scalar>&) {}
        virtual void operator/=(const faPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixedFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
