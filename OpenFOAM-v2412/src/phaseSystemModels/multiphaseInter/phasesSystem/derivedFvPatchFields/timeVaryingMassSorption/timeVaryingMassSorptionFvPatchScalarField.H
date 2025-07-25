/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::timeVaryingMassSorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a first order fixed-value condition
    for a given scalar field to model time-dependent adsorption-desoprtion
    processes to be used with the \c interfaceOxideRate mass model

    \f[
        \frac{d c}{d t} =
            k_{abs} w (c_{int} - c_{p_{w}}) + k_{des} (c_{p_{w}} - c_{int})
    \f]

    \f[
        w = \max(1 - c_{p_{w}}/max, 0)
    \f]

    where
    \vartable
      c_{int}      | Concentration at cell
      c_{p_{w}}    | Concentration at wall
      k_{abs}      | Adsorption rate constant [1/s]
      k_{des}      | Desorption rate constant [1/s]
      w            | Weight function
      max          | Max concentration at wall
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type        timeVaryingMassSorption;
        kbas        <scalar>;
        max         <scalar>;

        // Optional entries
        kdes        <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd   | Deflt
      type     | Type name: timeVaryingAdsorption | word   | yes    | -
      kbas     | Adsorption rate constant         | scalar | yes    | -
      max      | Maximum concentation at wall     | scalar | yes    | -
      kdes     | Desorption rate constant         | scalar | no     | 0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

SourceFiles
    timeVaryingMassSorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMassSorptionFvPatchScalarField_H
#define timeVaryingMassSorptionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class timeVaryingMassSorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMassSorptionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    // Public Enumeration

        //- Enumeration defining the available ddt schemes
        enum ddtSchemeType
        {
            tsEuler,
            tsCrankNicolson,
            tsBackward
        };


private:

    // Private Data

        //- Time scheme type names
        static const Enum<ddtSchemeType> ddtSchemeTypeNames_;

        //- Adsorption rate constant
        scalar kabs_;

        //- Maximum level of adsorption of a given substance on patch
        scalar max_;

        //- Desorption rate constant
        scalar kdes_;


public:

    //- Runtime type information
    TypeName("timeVaryingMassSorption");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMassSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMassSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- timeVaryingMassSorptionFvPatchScalarField onto a new patch
        timeVaryingMassSorptionFvPatchScalarField
        (
            const timeVaryingMassSorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMassSorptionFvPatchScalarField
        (
            const timeVaryingMassSorptionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        timeVaryingMassSorptionFvPatchScalarField
        (
            const timeVaryingMassSorptionFvPatchScalarField&,
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

        // Help

            //- Return source rate
            tmp<scalarField> source() const;


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
