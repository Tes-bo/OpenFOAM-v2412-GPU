/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::atmBoundaryLayerInletVelocityFvPatchVectorField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a log-law type ground-normal inlet
    boundary condition for the streamwise component of wind velocity
    (i.e. \c u) for homogeneous, two-dimensional, dry-air, equilibrium
    and neutral atmospheric boundary layer modelling.

    The ground-normal streamwise flow speed profile expression:

        \f[
            u = \frac{u^*}{\kappa} \ln \left( \frac{z - d + z_0}{z_0} \right)
        \f]

        \f[
            v = w = 0
        \f]

    where
    \vartable
      u        | Ground-normal streamwise flow speed profile    [m/s]
      v        | Spanwise flow speed                            [m/s]
      w        | Ground-normal flow speed                       [m/s]
      u^*      | Friction velocity                              [m/s]
      \kappa   | von Kármán constant                            [-]
      z        | Ground-normal coordinate component             [m]
      d        | Ground-normal displacement height              [m]
      z_0      | Aerodynamic roughness length                   [m]
    \endvartable

    Required fields:
    \verbatim
        U   | Velocity                                          [m/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        // Mandatory entries (unmodifiable)
        type            atmBoundaryLayerInletVelocity;

        // Mandatory/Optional (inherited) entries (unmodifiable)
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                              | Type | Reqd | Deflt
      type      | Type name: atmBoundaryLayerInletVelocity | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link atmBoundaryLayer.H \endlink
     - \link inletOutletFvPatchField.H \endlink

See also
    - Foam::atmBoundaryLayer
    - Foam::atmBoundaryLayerInletKFvPatchScalarField
    - Foam::atmBoundaryLayerInletEpsilonFvPatchScalarField
    - Foam::atmBoundaryLayerInletOmegaFvPatchScalarField

SourceFiles
    atmBoundaryLayerInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerInletVelocityFvPatchVectorField_H
#define atmBoundaryLayerInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmBoundaryLayerInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerInletVelocityFvPatchVectorField
:
    public inletOutletFvPatchVectorField,
    public atmBoundaryLayer
{

public:

    //- Runtime type information
    TypeName("atmBoundaryLayerInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmBoundaryLayerInletVelocityFvPatchVectorField onto a new patch
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const atmBoundaryLayerInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const atmBoundaryLayerInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
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
                const fvPatchVectorField&,
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
