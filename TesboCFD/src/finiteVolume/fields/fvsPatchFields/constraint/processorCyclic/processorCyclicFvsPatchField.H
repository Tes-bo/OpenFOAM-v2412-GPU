/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::processorCyclicFvsPatchField

Description
    This boundary condition enables processor communication across cyclic
    patches.

    The "value" entry is MUST_READ.

SourceFiles
    processorCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorCyclicFvsPatchField_H
#define Foam_processorCyclicFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "processorCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorCyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorCyclicFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const processorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(processorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary.
        //- The "value" entry is MUST_READ.
        processorCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        processorCyclicFvsPatchField
        (
            const processorCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        processorCyclicFvsPatchField(const processorCyclicFvsPatchField<Type>&);

        //- Construct as copy setting internal field reference
        processorCyclicFvsPatchField
        (
            const processorCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return fvsPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return fvsPatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~processorCyclicFvsPatchField() = default;


    // Member Functions

        //- True if running parallel
        virtual bool coupled() const
        {
            return UPstream::parRun();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
