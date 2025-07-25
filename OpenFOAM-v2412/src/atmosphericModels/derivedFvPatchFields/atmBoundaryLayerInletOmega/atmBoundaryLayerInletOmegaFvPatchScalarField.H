/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
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
    Foam::atmBoundaryLayerInletOmegaFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a log-law type ground-normal inlet
    boundary condition for the specific dissipation rate (i.e. \c omega)
    for homogeneous, two-dimensional, dry-air, equilibrium and neutral
    atmospheric boundary layer modelling.

    The ground-normal \c omega profile expression:

        \f[
            \omega = \frac{u^*}{\kappa \sqrt{C_\mu}} \frac{1}{z - d + z_0}
        \f]

    where
    \vartable
      \omega | Ground-normal specific dissipation rate profile [m^2/s^3]
      u^*    | Friction velocity                               [m/s]
      \kappa | von Kármán constant                             [-]
      C_\mu  | Empirical model constant                        [-]
      z      | Ground-normal coordinate component              [m]
      d      | Ground-normal displacement height               [m]
      z_0    | Aerodynamic roughness length                    [m]
    \endvartable

    Required fields:
    \verbatim
        omega   | Specific dissipation rate                    [1/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        // Mandatory entries (unmodifiable)
        type            atmBoundaryLayerInletOmega;

        // Mandatory/Optional (inherited) entries (unmodifiable)
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type      | Type name: atmBoundaryLayerInletOmega  | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link atmBoundaryLayer.H \endlink
     - \link inletOutletFvPatchField.H \endlink

See also
    - Foam::atmBoundaryLayer
    - Foam::atmBoundaryLayerInletVelocityFvPatchVelocityField
    - Foam::atmBoundaryLayerInletEpsilonFvPatchScalarField
    - Foam::atmBoundaryLayerInletKFvPatchScalarField

SourceFiles
    atmBoundaryLayerInletOmegaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerInletOmegaFvPatchScalarField_H
#define atmBoundaryLayerInletOmegaFvPatchScalarField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class atmBoundaryLayerInletOmegaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerInletOmegaFvPatchScalarField
:
    public inletOutletFvPatchScalarField,
    public atmBoundaryLayer
{

public:

    //- Runtime type information
    TypeName("atmBoundaryLayerInletOmega");


    // Constructors

        //- Construct from patch and internal field
        atmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmBoundaryLayerInletOmegaFvPatchScalarField onto a new patch
        atmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const atmBoundaryLayerInletOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        atmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const atmBoundaryLayerInletOmegaFvPatchScalarField&,
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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


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


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
