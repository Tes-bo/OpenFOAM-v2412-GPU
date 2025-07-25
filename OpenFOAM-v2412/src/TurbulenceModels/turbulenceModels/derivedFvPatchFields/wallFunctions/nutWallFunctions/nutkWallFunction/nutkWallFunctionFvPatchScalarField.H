/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::nutkWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the turbulent
    viscosity (i.e. \c nut) based on the turbulent kinetic energy,
    (i.e. \c k) for for low- and high-Reynolds number applications.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutkWallFunction;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type | Reqd | Deflt
      type      | Type name: nutkWallFunction   | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink
      - \link wallFunctionBlenders.H \endlink

SourceFiles
    nutkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkWallFunctionFvPatchScalarField_H
#define nutkWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFunctionBlenders.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField,
    private wallFunctionBlenders
{
protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutkWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutkWallFunctionFvPatchScalarField
        //- onto a new patch
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&,
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

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
