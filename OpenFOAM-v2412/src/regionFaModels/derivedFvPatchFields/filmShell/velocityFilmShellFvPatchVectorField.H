/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::compressible::velocityFilmShellFvPatchVectorField

Group
    grpLiquidFilmBoundaryConditions

Description
    This boundary condition provides a coupled velocity condition between
    a primary region (3D mesh) and a liquid-film model (2D mesh).

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                velocityFilmShell;
        liquidFilmModel     <word>;

        // Optional entries
        zeroWallVelocity    <bool>;

        // Inherited entries
        ...
        active              <bool>;
        infoOutput          <bool>;
        U                   <word>;
        pRef                <scalar>;
        T0                  <scalar>;
        deltaWet            <scalar>;
        h0                  <scalar>;

        thermo
        {
            H2O;
        }

        turbulence          laminar;
        laminarCoeffs
        {
            friction        ManningStrickler;   // Wall friction model
            n               0.005;              // Manning number
            Cf              0.9;                // Gas friction
        }

        injectionModels
        (
            curvatureSeparation
        );

        forces ();

        curvatureSeparationCoeffs
        {
            definedPatchRadii  0;
        }

        region              film;
        value               uniform (0 0 0);
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                        | Type | Reqd | Deflt
      type        | Type name: velocityFilmShell       | word | yes  | -
      liquidFilmModel     | Film model               | word   | yes  | -
      zeroWallVelocity    | Flag to fix zero U for primary flow <!--
                  -->                                  | bool | no   | true
      U           | Name of the primary U              | word | yes  | -
      pRef        | Reference pressure for thermo    | scalar | yes  | -
      T0          | Film initial temperature         | scalar | no   | READ
      thermo      | Flow thermo                     | wordRes | yes  | -
      turbulence  | Type of film turbulence model      | word | yes  | -
      injectionModels   | Lagrangian injection         |      | no   | -
      forces      | Film force models               | wordRes | no   | -
      deltaWet    | Wet film thickness               | scalar | no   | 1e-4
      h0          | Numerical minimum thickness      | scalar | no   | 1e-7
      region      | Name of the 2D region            | word   | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link liquidFilmBase.H \endlink
      - \link mixedFvPatchFields.H \endlink

SourceFiles
    velocityFilmShellFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef velocityFilmShellFvPatchVectorField_H
#define velocityFilmShellFvPatchVectorField_H

#include "liquidFilmBase.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class velocityFilmShellFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class velocityFilmShellFvPatchVectorField
:
    public mixedFvPatchField<vector>
{
    // Typedefs

        //- The finite-area region model
        typedef regionModels::areaSurfaceFilmModels::liquidFilmBase baffleType;


    // Private Data

        //- The liquid film model
        autoPtr<baffleType> baffle_;

        //- Dictionary
        mutable dictionary dict_;

        //- Time index to evolve the film
        label curTimeIndex_;

        //- Flag to set velocity to zero or film velocity
        bool zeroWallVelocity_;


public:

    //- Runtime type information
    TypeName("velocityFilmShell");


    // Constructors

        //- Construct from patch and internal field
        velocityFilmShellFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        velocityFilmShellFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- velocityFilmShellFvPatchVectorField onto a new patch
        velocityFilmShellFvPatchVectorField
        (
            const velocityFilmShellFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        velocityFilmShellFvPatchVectorField
        (
            const velocityFilmShellFvPatchVectorField&,
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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
