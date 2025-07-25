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
    Foam::speciesSorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a first-order zero-gradient
    condition for a given scalar field to model time-dependent
    adsorption-desorption processes.

    \f[
        \frac{d c}{d t} = k_{ads} (c_{eq} - c_{abs})
    \f]

    where
    \vartable
      c_{eq}       | Equilibrium concentration
      c_{abs}      | Absorbed at wall
      k_{ads}      | Adsorption rate constant [1/s]
    \endvartable

    \f[
        c_{eq} = c_{max} \frac{k_l \, c_{int}}{1 + k_l \, c_{int}}
    \f]

    where
    \vartable
      c_{max}      | Maximum concentration
      k_l          | Langmuir constant
      c_{int}      | Local cell value concentration
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                    speciesSorption;
        equilibriumModel        <word>;
        kinematicModel          <word>;
        kabs                    <scalar>;
        kl                      <scalar>;
        max                     <scalar>;
        thickness               <PatchFunction1<scalar>>;
        rhoS                    <scalar>;

        // Optional entries
        dfldp                   <scalarField>;
        mass                    <scalarField>;
        pName                   <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                        | Type   | Reqd | Deflt
      type     | Type name: speciesSorption         | word   | yes  | -
      equilibriumModel    | Equilibrium model       | word   | yes  | -
      kinematicModel      | Kinematic model         | word   | yes  | -
      kabs     | Adsorption rate constant [1/s]     | scalar | yes  | -
      kl       | Langmuir constant [1/Pa]           | scalar | yes  | -
      max      | Maximum concentation at wall [mol/kg] | scalar | yes | -
      thickness| Solid thickness along the patch   <!--
               -->                   | PatchFunction1\<scalar\> | yes | -
      rhoS     | Solid density                      | scalar | yes  | -
      dfldp    | Source on cells next to patch      | scalarField | no  | Zero
      mass     | Absorbed mass per kg of absorbent [mol/kg] <!--
               -->                                  | scalarField | no  | Zero
      pName    | Name of operand pressure field     | word   | no   | p
    \endtable

    Options for the \c equilibriumModel entry:
    \verbatim
      Langmuir    | Langmuir model
    \endverbatim

    Options for the \c kinematicModel entry:
    \verbatim
      PseudoFirstOrder    | Pseudo first-order model
    \endverbatim

    The inherited entries are elaborated in:
      - \link zeroGradientFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    speciesSorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_speciesSorptionFvPatchScalarField_H
#define Foam_speciesSorptionFvPatchScalarField_H

#include "boundarySourcePatch.H"
#include "zeroGradientFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class speciesSorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class speciesSorptionFvPatchScalarField
:
    public zeroGradientFvPatchField<scalar>,
    public boundarySourcePatch
{
public:

    // Public Enumeration

        //- Options for the equilibrum model
        enum equilibriumModelType : char
        {
            LANGMUIR = 0
        };

        //- Options for the kinematic model
        enum kineticModelType : char
        {
            PseudoFirstOrder = 0
        };

        //- Names for equilibriumModelType
        static const Enum<equilibriumModelType> equilibriumModelTypeNames;

        //- Names for kineticModelType
        static const Enum<kineticModelType> kinematicModelTypeNames;


private:

    // Private Data

        //- Equilibrium model
        enum equilibriumModelType equilibriumModel_;

        //- Kinematic model
        enum kineticModelType kinematicModel_;

        //- Solid thickness along the patch
        autoPtr<PatchFunction1<scalar>> thicknessPtr_;

        //- Adsorption rate constant  [1/sec]
        scalar kabs_;

        //- Langmuir adsorption constant [1/Pa]
        scalar kl_;

        //- Maximum density on patch [mol/kg]
        scalar max_;

        //- Solid density
        scalar rhoS_;

        //- Name of operand pressure field
        word pName_;

        //- Source on cells next to patch [mol/kg/sec]
        scalarField dfldp_;

        //- Absorbed mass per kg of absorbent [mol/kg]
        scalarField mass_;


    // Private Member Functions

        //- Calculate the mole fraction fields
        tmp<scalarField> calcMoleFractions() const;

        //- Lookup (or create) field for output
        volScalarField& field(const word&, const dimensionSet&) const;


public:

    //- Runtime type information
    TypeName("speciesSorption");


    // Constructors

        //- Construct from patch and internal field
        speciesSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        speciesSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- speciesSorptionFvPatchScalarField onto a new patch
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&,
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

        //- Access to mass
        tmp<scalarField> mass() const;

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
