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
    Foam::atmAlphatkWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the kinematic
    turbulent thermal conductivity (i.e. \c alphat) for atmospheric boundary
    layer modelling.  It assumes a logarithmic distribution of the potential
    temperature within the first cell.

    Required fields:
    \verbatim
      alphat    | Kinematic turbulent thermal conductivity    [m2/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmAlphatkWallFunction;
        Pr              <Function1<scalar>>;
        Prt             <PatchFunction1<scalar>>;
        z0              <PatchFunction1<scalar>>;

        // Optional entries
        Cmu             <scalar>;
        kappa           <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                       | Type   | Reqd     | Deflt
      type     | Type name: atmAlphatkWallFunction | word   | yes      | -
      Pr       | Molecular Prandtl number     | Function1<scalar>  | yes | -
      Prt      | Turbulent Prandtl number     | PatchFunction1<scalar> | yes | -
      z0       | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
      Cmu      | Empirical model constant         | scalar | no        | 0.09
      kappa    | von Kármán constant              | scalar | no        | 0.41
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink
      - \link Function1.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    atmAlphatkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmAlphatkWallFunctionFvPatchScalarField_H
#define atmAlphatkWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
       Class atmAlphatkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmAlphatkWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected Data

        //- Empirical model constant
        const scalar Cmu_;

        //- von Kármán constant
        const scalar kappa_;

        //- Molecular Prandtl number
        autoPtr<Function1<scalar>> Pr_;

        //- Turbulent Prandtl number field
        autoPtr<PatchFunction1<scalar>> Prt_;

        //- Surface roughness length [m]
        autoPtr<PatchFunction1<scalar>> z0_;


        //- Solution parameters

            static scalar tolerance_;
            static label maxIters_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("atmAlphatkWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmAlphatkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmAlphatkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmAlphatkWallFunctionFvPatchScalarField
        //- onto a new patch
        atmAlphatkWallFunctionFvPatchScalarField
        (
            const atmAlphatkWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmAlphatkWallFunctionFvPatchScalarField
        (
            const atmAlphatkWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmAlphatkWallFunctionFvPatchScalarField
        (
            const atmAlphatkWallFunctionFvPatchScalarField&,
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
