/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::nutWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    The class \c nutWallFunction is an abstract base class that
    hosts \f$y^+\f$ calculation methods and common functions for
    \c nut wall-function boundary conditions.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Top-level entries
        ...

        // Optional entries
        U             U;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Reqd  | Deflt
      U         | Name of operand velocity field | word   | no    | U
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink
      - \link wallFunctionCoefficients.H \endlink

SourceFiles
    nutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutWallFunctionFvPatchScalarField_H
#define nutWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "wallFunctionCoefficients.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class nutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of velocity field
        //  Default is null (not specified) in which case the velocity is
        //  retrieved from the turbulence model
        word UName_;

        //- Wall-function coefficients
        wallFunctionCoefficients wallCoeffs_;


    // Protected Member Functions

        //- Helper to return the velocity field either from the turbulence
        //- model (default) or the mesh database
        virtual const volVectorField& U(const turbulenceModel& turb) const;

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const = 0;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutWallFunctionFvPatchScalarField
        //- onto a new patch
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        // No clone methods - abstract class


    // Member Functions

        //- Return the nut patchField for the given wall patch
        static const nutWallFunctionFvPatchScalarField& nutw
        (
            const turbulenceModel& turbModel,
            const label patchi
        );

        //- Calculate and return the yPlus at the boundary
        //  yPlus is the first-cell-centre height from boundary in wall units
        virtual tmp<scalarField> yPlus() const = 0;

        //- Return wallFunctionCoefficients
        inline const wallFunctionCoefficients& wallCoeffs() const noexcept
        {
            return wallCoeffs_;
        }


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
