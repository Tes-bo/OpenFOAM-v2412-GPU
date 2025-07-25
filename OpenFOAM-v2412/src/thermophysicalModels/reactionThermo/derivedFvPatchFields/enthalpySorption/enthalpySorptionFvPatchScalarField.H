/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::enthalpySorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This is a temperature boundary condition which works
    in conjunction with the \c speciesSorption condition for species.

    This boundary condition substracts or adds enthalpy associated with the
    adsorption provided by the \c speciesSorption condition.

    It can handle two enthalpy models:

    1) Estimate
    2) Calculated

    On top of this, the enthalpy associated with the sensible enthalpy
    corresponding with the species transport can be added using \c includeHs.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                enthalpySorption;
        enthalpyModel       <word>;
        species             <word>;

        // Conditional mandatory entries

        // when enthalpyModel == calculated
        enthalpyTable       <Function1<scalar>>

            // enthalpyTable
            // {
            //     type            table;
            //     values          ((0 0)(1 50));
            // }

        // Optional entries
        includeHs           <bool>;
        C                   <scalar>;
        Hvap                <scalar>;
        dhdt                <scalarField>;
        p                   <word>;
        T                   <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                        | Type   | Reqd   | Deflt
      type     | Type name: enthalpySorption        | word   | yes    | -
      enthalpyModel | Adsorption enthalpy model     | word   | yes    | -
      species  | Name of associated species         | word   | yes    | -
      enthalpyTable | Calculated enthalpy model table       <!--
               -->                         | Function1\<scalar\> | no | -
      includeHs | Include sensible enthalpy         | bool   | no     | true
      C        | Estimate enthalpy model constant   | scalar | no     | 0
      Hvap     | Evaporation heat for species       | scalar | no     | 0
      p        | Name of operand pressure field     | word   | no     | p
      T        | Name of operand temperature field  | word   | no     | T
      dhdt     | Enthalpy change on cells next to patch [J/kg]       <!--
               -->                             | scalarField | no     | Zero
    \endtable

    Options for the \c enthalpyModel entry:
    \verbatim
      estimated    | Enthalpy is estimated
      calculated   | Enthalpy is calculated based on enthalpyTable
    \endverbatim

    The inherited entries are elaborated in:
      - \link zeroGradientFvPatchFields.H \endlink
      - \link Function1.H \endlink

SourceFiles
    enthalpySorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_enthalpySorptionFvPatchScalarField_H
#define Foam_enthalpySorptionFvPatchScalarField_H

#include "boundarySourcePatch.H"
#include "zeroGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class enthalpySorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class enthalpySorptionFvPatchScalarField
:
    public zeroGradientFvPatchField<scalar>,
    public boundarySourcePatch
{
    // Private Enumeration

        //- Options for the enthalpy model
        enum enthalpyModelType : char
        {
            estimated,
            calculated
        };

        //- Names for enthalpyModelType
        static const Enum<enthalpyModelType> enthalpyModelTypeNames;


    // Private Data

        //- Enthalpy model
        enum enthalpyModelType enthalpyModel_;

        //- Include sensible enthalpy of the species
        bool includeHs_;

        //- Load enthalpy table for calculated model
        autoPtr<Function1<scalar>> enthalpyMassLoadPtr_;

        //- Estimated enthalpy model constant
        scalar C_;

        //- Heat of evaporation of species
        scalar Hvap_;

        //- Name of operand species field
        word speciesName_;

        //- Name of operand pressure field
        word pName_;

        //- Name of operand temperature field
        word TName_;

        //- Enthalpy change on cells next to patch [J/kg]
        scalarField dhdt_;


public:

    //- Runtime type information
    TypeName("enthalpySorption");


    // Constructors

        //- Construct from patch and internal field
        enthalpySorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        enthalpySorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- enthalpySorptionFvPatchScalarField onto a new patch
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
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
            const fvPatchScalarField&,
            const labelList&
        );


    // Evaluation

        //- Source of cells next to the patch
        virtual tmp<scalarField> patchSource() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


    // I-O

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
