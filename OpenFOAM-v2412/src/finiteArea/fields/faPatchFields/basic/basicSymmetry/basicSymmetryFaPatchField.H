/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::basicSymmetryFaPatchField

Description
    A symmetry patch

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    basicSymmetryFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_basicSymmetryFaPatchField_H
#define Foam_basicSymmetryFaPatchField_H

#include "transformFaPatchField.H"
#include "symmetryFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class basicSymmetryFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class basicSymmetryFaPatchField
:
    public transformFaPatchField<Type>
{
public:

    // Constructors

        //- Construct from patch and internal field
        basicSymmetryFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        //  The "value" entry: NO_READ
        basicSymmetryFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given basicSymmetryFaPatchField
        //- onto a new patch
        basicSymmetryFaPatchField
        (
            const basicSymmetryFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Default copy construct
        basicSymmetryFaPatchField
        (
            const basicSymmetryFaPatchField<Type>&
        ) = default;

        //- Construct as copy setting internal field reference
        basicSymmetryFaPatchField
        (
            const basicSymmetryFaPatchField<Type>&,
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

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            //  Default argument needed to allow call in constructors
            //  HJ, 30/Jun/2009
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
            );

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const;


    // Member Operators

        //- Inherit assignment operators
        using transformFaPatchField<Type>::operator=;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> basicSymmetryFaPatchField<scalar>::snGrad() const;

template<>
void basicSymmetryFaPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicSymmetryFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
