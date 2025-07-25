/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::compressible::
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField

Description
    This boundary condition provides a thermal wall function for turbulent
    thermal diffusivity (usually\c alphat) based on the Jayatilleke model for
    the Eulerian multiphase solvers.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            alphatPhaseChangeJayatillekeWallFunction;

        // Optional entries
        Prt             <scalar>;
        Cmu             <scalar>;
        kappa           <scalar>;
        E               <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      type      | Type name:                             <!--
      -->  compressible::alphatPhaseChangeJayatillekeWallFunction <!--
      -->                                                 | word | yes  | -
      Prt       | Turbulent Prandtl number              | scalar | no   | 0.85
      Cmu       | Empirical model coefficient           | scalar | no   | 0.09
      kappa     | Von Karman constant                   | scalar | no   | 0.41
      E         | Wall roughness parameter              | scalar | no   | 9.8
    \endtable

    The inherited entries are elaborated in:
      - \link alphatPhaseChangeWallFunctionFvPatchScalarField.H \endlink

See also
    Foam::compressible::alphatPhaseChangeWallFunctionFvPatchScalarField

SourceFiles
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField_H
#define compressible_alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
 Class alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Empirical model coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- Wall roughness parameter
        scalar E_;

        // Solution parameters

            //- Absolute tolerance
            static scalar tolerance_;

            //- Maximum number of iterations
            static label maxIters_;


    // Protected Member Functions

        //- Check the type of the patch
        void checkType();

        //- 'P' function
        tmp<scalarField> Psmooth(const scalarField& Prat) const;

        //- Calculate y+ at the edge of the thermal laminar sublayer
        tmp<scalarField> yPlusTherm
        (
            const scalarField& P,
            const scalarField& Prat
        ) const;

        //- Update turbulent thermal diffusivity
        tmp<scalarField> calcAlphat
        (
            const scalarField& prevAlphat
        ) const;


public:

    //- Runtime type information
    TypeName("compressible::alphatPhaseChangeJayatillekeWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        //- onto a new patch
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField&,
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

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
