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
    Foam::calculatedFaePatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    calculatedFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_calculatedFaePatchField_H
#define Foam_calculatedFaePatchField_H

#include "faePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class calculatedFaePatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class calculatedFaePatchField
:
    public faePatchField<Type>
{
public:

    //- Runtime type information
    TypeName("calculated");  // fieldTypes::calculatedTypeName_()


    // Constructors

        //- Construct from patch and internal field
        calculatedFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        calculatedFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary& dict,
            //! The "value" entry (default: mandatory)
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Construct by mapping given patch field onto a new patch
        calculatedFaePatchField
        (
            const calculatedFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        calculatedFaePatchField(const calculatedFaePatchField<Type>&);

        //- Construct as copy setting internal field reference
        calculatedFaePatchField
        (
            const calculatedFaePatchField<Type>&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Return clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return faePatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, edgeMesh>& iF
        ) const
        {
            return faePatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~calculatedFaePatchField() = default;


    // Member Functions

        //- True: the patch field fixes a value.
        virtual bool fixesValue() const { return true; }


    // Other

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "calculatedFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
