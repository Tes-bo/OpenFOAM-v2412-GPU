/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::nutUBlendedWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the turbulent
    viscosity (i.e. \c nut) based on velocity (i.e. \c U) using a
    binomial-function wall-function blending method between the viscous
    and inertial sublayer predictions of \c nut for low- and high-Reynolds
    number applications.

    \f[
        u_\tau = (u_{\tau,vis}^n + u_{\tau,log}^n)^{1/n}
    \f]

    where
    \vartable
      u_\tau     | Friction velocity
      u_{\tau,vis} | Friction velocity in the viscous sublayer
      u_{\tau,log} | Friction velocity in the inertial sublayer
    \endvartable

    Reference:
    \verbatim
        See the section that describes 'automatic wall treatment':
            Menter, F., Ferreira, J. C., Esch, T., Konno, B. (2003).
            The SST turbulence model with improved wall treatment
            for heat transfer predictions in gas turbines.
            In Proceedings of the International Gas Turbine Congress.
            November, 2003. Tokyo, Japan. pp. 2-7.
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutUBlendedWallFunction;

        // Optional entries
        n               4.0;

        // Inherited entries
        ...
     }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type   | Reqd | Deflt
      type      | Type name: nutUBlendedWallFunction  | word   | yes  | -
      n         | Blending factor                     | scalar | no   | 4.0
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink

Note
    - The full 'automatic wall treatment' description also requires use of the
    \link omegaWallFunctionFvPatchScalarField.H \endlink with the \c blending
    option \c binomial or with the deprecated \c blended flag set to \c on.
    - Suffers from non-exact restart since \c correctNut() (called through
    \c turbulence->validate) returns a slightly different value every time
    it is called.
    See \link nutUSpaldingWallFunctionFvPatchScalarField.C \endlink.

SourceFiles
    nutUBlendedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUBlendedWallFunctionFvPatchScalarField_H
#define nutUBlendedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUBlendedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUBlendedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Model coefficient; default = 4
        scalar n_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutUBlendedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutUBlendedWallFunctionFvPatchScalarField
        //- onto a new patch
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&,
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

        // Evaluation

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
