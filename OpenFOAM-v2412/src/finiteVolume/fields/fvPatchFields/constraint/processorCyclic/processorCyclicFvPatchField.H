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
    Foam::processorCyclicFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enables processor communication across cyclic
    patches.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            processor;
    }
    \endverbatim

See also
    Foam::processorFvPatchField

SourceFiles
    processorCyclicFvPatchField.C
    processorCyclicFvPatchFields.H
    processorCyclicFvPatchFields.C
    processorCyclicFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorCyclicFvPatchField_H
#define Foam_processorCyclicFvPatchField_H

#include "processorCyclicFvPatch.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class processorCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorCyclicFvPatchField
:
    public processorFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const processorCyclicFvPatch& procPatch_;

public:

    //- Runtime type information
    TypeName(processorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorCyclicFvPatchField onto a
        //  new patch
        processorCyclicFvPatchField
        (
            const processorCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        processorCyclicFvPatchField(const processorCyclicFvPatchField<Type>&);

        //- Construct as copy setting internal field reference
        processorCyclicFvPatchField
        (
            const processorCyclicFvPatchField<Type>&,
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
    virtual ~processorCyclicFvPatchField() = default;


    // Member Functions

        // Access

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
