/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::compressible::RASModels::alphatFilmWallFunctionFvPatchScalarField

Group
    grpSurfaceFilmBoundaryConditions grpCmpWallFunctions

Description
    This boundary condition provides a turbulent thermal diffusivity condition
    when using wall functions, for use with surface film models.  This
    condition varies from the standard wall function by taking into account any
    mass released from the film model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            alphatFilmWallFunction;

        // Optional entries
        filmRegion      <word>;
        B               <scalar>;
        yPlusCrit       <scalar>;
        Cmu             <scalar>;
        kappa           <scalar>;
        Prt             <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                    | Type   | Reqd  | Deflt
      type       | Type name: alphatFilmWallFunction | word | yes  | -
      filmRegion | Name of film region   | word | no | surfaceFilmProperties
      B          | Model coefficient              | scalar | no    | 5.5
      yPlusCrit  | Critical y+ for transition to turbulence | scalar | no | 11.5
      Cmu        | Empirical model coefficient    | scalar | no    | 0.09
      Prt        | Turbulent Prandtl number       | scalar | no    | 0.85
      kappa      | von Karman constant            | scalar | no    | 0.41
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

SourceFiles
    alphatFilmWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatFilmWallFunctionFvPatchScalarField_H
#define alphatFilmWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
          Class alphatFilmWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatFilmWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of film region
        word filmRegionName_;

        //- B Coefficient
        scalar B_;

        //- y+ value for laminar -> turbulent transition
        scalar yPlusCrit_;

        //- Empirical model coefficient
        scalar Cmu_;

        //- Von-Karman constant
        scalar kappa_;

        //- Turbulent Prandtl number
        scalar Prt_;


public:

    //- Runtime type information
    TypeName("alphatFilmWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- alphatFilmWallFunctionFvPatchScalarField
        //- onto a new patch
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        alphatFilmWallFunctionFvPatchScalarField
        (
            const alphatFilmWallFunctionFvPatchScalarField&,
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

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
