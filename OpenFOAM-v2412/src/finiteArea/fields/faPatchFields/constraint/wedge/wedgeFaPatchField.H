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
    Foam::wedgeFaPatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    wedgeFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wedgeFaPatchField_H
#define Foam_wedgeFaPatchField_H

#include "transformFaPatchField.H"
#include "wedgeFaPatch.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class wedgeFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wedgeFaPatchField
:
    public transformFaPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(wedgeFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        wedgeFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        //  The "value" entry: NO_READ
        wedgeFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wedgeFaPatchField onto a new patch
        wedgeFaPatchField
        (
            const wedgeFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        wedgeFaPatchField
        (
            const wedgeFaPatchField<Type>&,
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

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
            );


            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> wedgeFaPatchField<scalar>::snGrad() const;

template<>
void wedgeFaPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wedgeFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
