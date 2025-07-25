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
    Foam::processorFaePatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    processorFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorFaePatchField_H
#define Foam_processorFaePatchField_H

#include "coupledFaePatchField.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class processorFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFaePatchField
:
    public coupledFaePatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const processorFaPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(processorFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorFaePatchField onto a new patch
        processorFaePatchField
        (
            const processorFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, edgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        processorFaePatchField(const processorFaePatchField<Type>&);

        //- Construct as copy setting internal field reference
        processorFaePatchField
        (
            const processorFaePatchField<Type>&,
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
    virtual ~processorFaePatchField() = default;


    // Member Functions

        //- Return true if running parallel
        virtual bool coupled() const
        {
            return UPstream::parRun();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
