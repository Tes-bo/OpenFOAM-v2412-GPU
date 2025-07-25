/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::compressible::RASModels::nutkFilmWallFunctionFvPatchScalarField

Group
    grpSurfaceFilmBoundaryConditions grpCmpWallFunctions

Description
    This boundary condition provides a turbulent viscosity condition when
    using wall functions, based on turbulence kinetic energy, for use with
    surface film models.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutkFilmWallFunction;

        // Optional entries
        filmRegion      <word>;
        B               <scalar>;
        yPlusCrit       <scalar>;

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
    \endtable

    The inherited entries are elaborated in:
      - \link nutkWallFunctionFvPatchScalarField.H \endlink

SourceFiles
    nutkFilmWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkFilmWallFunctionFvPatchScalarField_H
#define nutkFilmWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
         Class nutkFilmWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkFilmWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of film region
        word filmRegionName_;

        //- B Coefficient
        scalar B_;

        //- y+ value for laminar -> turbulent transition
        scalar yPlusCrit_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutkFilmWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkFilmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutkFilmWallFunctionFvPatchScalarField
        //- onto a new patch
        nutkFilmWallFunctionFvPatchScalarField
        (
            const nutkFilmWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutkFilmWallFunctionFvPatchScalarField
        (
            const nutkFilmWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutkFilmWallFunctionFvPatchScalarField
        (
            const nutkFilmWallFunctionFvPatchScalarField&,
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

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
