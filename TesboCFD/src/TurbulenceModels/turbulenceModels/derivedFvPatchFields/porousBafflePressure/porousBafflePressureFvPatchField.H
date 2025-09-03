/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::porousBafflePressureFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition,
    using the \c cyclic condition as a base.

    The porous baffle introduces a pressure jump defined by:

    \f[
        \Delta p = -(D \mu U + 0.5 I \rho |U|^2 )L
    \f]

    where
    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m3]
        \mu    | laminar viscosity [Pa s]
        D      | Darcy coefficient
        I      | inertial coefficient
        L      | porous media length in the flow direction
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            porousBafflePressure;
        D               0.001;
        I               1000000;
        length          0.1;

        // Optional entries
        phi             phi;
        rho             rho;
        uniformJump     false;
        value           uniform 0;

        // Inherited entries
        patchType       cyclic;
        jump            uniform 0;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                          | Type | Reqd | Deflt
      type      | Type name: porousBafflePressure      | word | yes  | -
      D         | Darcy coefficient       | Function1<scalar> | yes  | -
      I         | Inertial coefficient    | Function1<scalar> | yes  | -
      length    | Porous media length in the flow direction   <!--
                -->                                 | scalar  | yes  | -
      uniformJump | Flag to apply a uniform pressure drop on  <!--
                --> the patch based on the             | bool | no   | false
      phi       | Name of flux field                   | word | no   | phi
      rho       | Name of density field                | word | no   | rho
      patchType | Underlying patch type should be \c cyclic | word | yes | -
      jump      | Jump value                    | scalarField | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link fixedJumpFvPatchField.H \endlink
      - \link Function1.H \endlink

Note
    The underlying \c patchType should be set to \c cyclic.

SourceFiles
    porousBafflePressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureFvPatchField_H
#define porousBafflePressureFvPatchField_H

#include "fixedJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class porousBafflePressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class porousBafflePressureFvPatchField
:
    public fixedJumpFvPatchField<scalar>
{
    // Private Data

        //- Name of flux field (default = phi)
        const word phiName_;

        //- Name of density field (default = rho)
        const word rhoName_;

        //- Darcy pressure loss coefficient
        autoPtr<Function1<scalar>> D_;

        //- Inertia pressure lost coefficient
        autoPtr<Function1<scalar>> I_;

        //- Porous media length
        scalar length_;

        //- Aplies uniform pressure drop
        bool uniformJump_;


public:

    //- Runtime type information
    TypeName("porousBafflePressure");


    // Constructors

        //- Construct from patch and internal field
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool needValue = true
        );

        //- Construct by mapping given porousBafflePressureFvPatchField
        //- onto a new patch
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&
        );

        //- Construct as copy setting internal field reference
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField&,
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
