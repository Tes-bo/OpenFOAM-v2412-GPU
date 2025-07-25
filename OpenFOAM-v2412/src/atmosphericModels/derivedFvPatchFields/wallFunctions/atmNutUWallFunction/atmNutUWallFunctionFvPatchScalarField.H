/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
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
    Foam::atmNutUWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    viscosity (i.e. \c nut) based on velocity (i.e. \c U) for atmospheric
    boundary layer modelling.  It is designed to be used in conjunction
    with the \c atmBoundaryLayerInletVelocity boundary condition.

    The governing equation of the boundary condition:

        \f[
            u = \frac{u^*}{\kappa} ln \left(\frac{z + z_0}{z_0}\right)
        \f]

    where
    \vartable
        u^*      | Friction velocity
        \kappa   | von Kármán constant
        z_0      | Surface roughness length [m]
        z        | Ground-normal coordinate
    \endvartable

    Required fields:
    \verbatim
      nut     | Turbulent viscosity         [m2/s]
      U       | Velocity                    [m/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmNutUWallFunction;
        z0              <PatchFunction1<scalar>>;

        // Optional entries
        boundNut        true;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd  | Deflt
      type     | Type name: atmNutUWallFunction   | word   | yes   | -
      z0       | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
      boundNut | Flag: zero-bound nut near wall   | bool   | no    | true
    \endtable

    The inherited entries are elaborated in:
      - \link nutUWallFunctionFvPatchScalarField.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    atmNutUWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmNutUWallFunctionFvPatchScalarField_H
#define atmNutUWallFunctionFvPatchScalarField_H

#include "nutUWallFunctionFvPatchScalarField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class atmNutUWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmNutUWallFunctionFvPatchScalarField
:
    public nutUWallFunctionFvPatchScalarField
{
    // Private Data

        //- Flag to zero-bound nut to prevent negative nut
        //- at the wall arising from negative heat fluxes
        const bool boundNut_;

        //- Surface roughness length field [m]
        autoPtr<PatchFunction1<scalar>> z0_;


protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("atmNutUWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmNutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmNutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmNutUWallFunctionFvPatchScalarField
        //- onto a new patch
        atmNutUWallFunctionFvPatchScalarField
        (
            const atmNutUWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmNutUWallFunctionFvPatchScalarField
        (
            const atmNutUWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmNutUWallFunctionFvPatchScalarField
        (
            const atmNutUWallFunctionFvPatchScalarField&,
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
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
