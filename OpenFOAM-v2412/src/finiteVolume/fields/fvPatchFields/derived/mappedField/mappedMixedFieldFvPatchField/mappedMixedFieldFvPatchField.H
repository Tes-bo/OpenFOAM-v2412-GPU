/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::mappedMixedFieldFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides a self-contained version of e.g. mapped
    boundary conditions
    - it expects two-way coupling (so the sampled field needs to run the
      same bc)
    - it does not use information on the patch; instead it holds the coupling
      data locally.

Usage
    \table
        Property     | Description                      | Required | Default
        field        | Name of field to be mapped       | no  | this field name
        weightField  | Name of field to use as weight   | no |
        setAverage   | Use average value                | no  | false |
        average      | Average value to use if \c setAverage = yes | partly |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedMixedField;
        field           T;              // optional field name
        setAverage      no;             // apply an average value
        average         0;              // average to apply if setAverage
        value           uniform 0;      // place holder
    }
    \endverbatim

Note
    Supports multi-world operation:

    +----+                      +----+
    |    |                      |    |
    |    |                      |    |
    +----+                      +----+
    worldA                      worldB
    regionA                     regionB
    patchA                      patchB


See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::mappedFieldFvPatchField
    Foam::mixedFvPatchField

SourceFiles
    mappedMixedFieldFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMixedFieldFvPatchField_H
#define mappedMixedFieldFvPatchField_H

#include "mappedPatchBase.H"
#include "mappedPatchFieldBase.H"
#include "mixedFvPatchFields.H"
#include "interpolation.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedMixedFieldFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedMixedFieldFvPatchField
:
    public mixedFvPatchField<Type>,
    public mappedPatchBase,
    public mappedPatchFieldBase<Type>
{

        //- Name of weight field to sample
        word weightFieldName_;


public:

    //- Runtime type information
    TypeName("mappedMixedField");


    // Constructors

        //- Construct from patch and internal field
        mappedMixedFieldFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedMixedFieldFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedMixedFieldFvPatchField
        //  onto a new patch
        mappedMixedFieldFvPatchField
        (
            const mappedMixedFieldFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedMixedFieldFvPatchField
        (
            const mappedMixedFieldFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        mappedMixedFieldFvPatchField
        (
            const mappedMixedFieldFvPatchField<Type>&,
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


    // Member functions

        // Mapping functions

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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedMixedFieldFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
