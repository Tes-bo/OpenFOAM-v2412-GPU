/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::outletMappedUniformInletFvPatchField

Group
    grpInletBoundaryConditions

Description
    The \c outletMappedUniformInlet is an inlet boundary condition that
    - averages patch fields of specified "outlet" patches
    and uniformly applies the averaged value over a specified inlet patch.
    - optionally, the averaged value can be scaled
    and/or offset by a specified value,
    and/or mapped by a specified time delay.

    The governing equation of the boundary condition is:

    \f[
        \phi_{inlet} =
            \sum_i \left( f_i \phi_{{outlet}_i} + \phi_{{offset}_i} \right)
    \f]

    where
    \vartable
      \phi_{inlet}   | Spatially-uniform patch-field value at an inlet patch
      \phi_{outlet}  | Averaged patch-field value at an outlet patch
      f              | User-defined fraction value
      \phi_{offset}  | User-defined offset value
      i              | Outlet-patch index
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            outletMappedUniformInlet;

        outlets
        {
            <outletName.1>
            {
                fraction    <Function1<scalar>>;
                offset      <Function1<Type>>;
                timeDelay   <Function1<scalar>>;
            }
            <outletName.2>
            {
                fraction    <Function1<scalar>>;
                offset      <Function1<Type>>;
                timeDelay   <Function1<scalar>>;
            }
            ...
        }

        // Optional entries
        uniformValue    <PatchFunction1<Type>>;
        phi             phi;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                          | Type   | Reqd | Deflt
      type        | Type name: outletMappedUniformInlet  | word   | yes  | -
      outlets     | Dictionary name: outlets             | dict   | yes  | -
      fraction    | Fraction value          | Function1\<scalar\> | no   | 1
      offset      | Offset value            | Function1\<Type\>   | no   | Zero
      timeDelay   | Time delay              | Function1\<scalar\> | no   | 0
      uniformValue | Base inlet patch field | PatchFunction1\<Type\> | no | Zero
      phi         | Name of operand flux field           | word   | no   | phi
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink
      - \link Function1.H \endlink

Note
  - Any negative input of \c timeDelay entry
    is forced to be zero without emitting any warnings.

See also
  - Foam::fixedValueFvPatchField
  - Foam::outletMappedUniformInletHeatAdditionFvPatchField
  - Foam::outletMappedUniformInletTemperatureFvPatchField

SourceFiles
    outletMappedUniformInletFvPatchField.C
    outletMappedUniformInletFvPatchFields.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletFvPatchField_H
#define outletMappedUniformInletFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class outletMappedUniformInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletMappedUniformInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Base inlet patch field
        autoPtr<PatchFunction1<Type>> uniformValuePtr_;

        //- List of outlet-patch names
        wordList outletNames_;

        //- List of outlet-patch field offsets
        PtrList<Function1<Type>> offsets_;

        //- List of outlet-patch field fractions
        PtrList<Function1<scalar>> fractions_;

        //- List of outlet-patch field time delays
        PtrList<Function1<scalar>> timeDelays_;

        //- List of outlet-patch mapping fields
        List<DynamicList<Type>> mapFields_;

        //- List of outlet-patch mapping times
        List<DynamicList<scalar>> mapTimes_;

        //- Name of operand flux field
        word phiName_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInlet");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletMappedUniformInletFvPatchField
        //- onto a new patch
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
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

        // Mapping

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

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "outletMappedUniformInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
