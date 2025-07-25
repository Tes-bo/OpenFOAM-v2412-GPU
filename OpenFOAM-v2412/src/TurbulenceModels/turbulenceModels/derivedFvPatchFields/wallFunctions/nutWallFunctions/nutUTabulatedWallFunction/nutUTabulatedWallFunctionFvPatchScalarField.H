/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::nutUTabulatedWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    viscosity (i.e. \c nut) based on velocity (i.e. \c U), for low- and
    high-Reynolds number applications.

    As input, the user specifies a look-up table of \c u+ as a function of
    near-wall Reynolds number.

    The table should be located in the \c $FOAM_CASE/constant directory.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutTabulatedWallFunction;
        uPlusTable      myUPlusTable;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                          | Type | Reqd | Deflt
      type        | Type name: nutUTabulatedWallFunction | word | yes  | -
      uPlusTable  | u+ as a function of Re table name    | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink
      - \link uniformInterpolationTable.H \endlink

Note
  - The tables are not registered since the same table object may be used for
    more than one patch.

SourceFiles
    nutUTabulatedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUTabulatedWallFunctionFvPatchScalarField_H
#define nutUTabulatedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "uniformInterpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class nutUTabulatedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUTabulatedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of u+ table
        word uPlusTableName_;

        //- u+ table
        uniformInterpolationTable<scalar> uPlusTable_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate wall u+ from table
        tmp<scalarField> calcUPlus(const scalarField& Rey) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutTabulatedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutUTabulatedWallFunctionFvPatchScalarField
        //- onto a new patch
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutUTabulatedWallFunctionFvPatchScalarField
        (
            const nutUTabulatedWallFunctionFvPatchScalarField&,
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
