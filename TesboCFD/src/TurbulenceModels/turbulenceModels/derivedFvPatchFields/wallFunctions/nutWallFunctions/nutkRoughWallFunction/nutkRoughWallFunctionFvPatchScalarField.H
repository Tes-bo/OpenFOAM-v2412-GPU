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
    Foam::nutkRoughWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall-function for the turbulent
    viscosity (i.e. \c nut) when using wall functions for rough walls,
    based on the turbulent kinetic energy (i.e. \c k). The condition
    manipulates the wall roughness parameter (i.e. \c E) to account
    for roughness effects.

    Parameter ranges:
    - roughness height (i.e. \c Ks)= sand-grain roughness (0 for smooth walls)
    - roughness constant (i.e. \c Cs) = 0.5-1.0

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutkRoughWallFunction;
        Ks              uniform 0;
        Cs              uniform 0.5;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type   | Reqd  | Dflt
      type        | Type name: nutkRoughWallFunction  | word   | yes   | -
      Ks          | Sand-grain roughness height  | scalarField | yes   | -
      Cs          | Roughness constant           | scalarField | yes   | -
    \endtable

    The inherited entries are elaborated in:
      - \link nutkWallFunctionFvPatchScalarField.H \endlink

SourceFiles
    nutkRoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkRoughWallFunctionFvPatchScalarField_H
#define nutkRoughWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nutkRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkRoughWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Roughness height
        scalarField Ks_;

        //- Roughness constant
        scalarField Cs_;


    // Protected Member Functions

        //- Compute the roughness function
        virtual scalar fnRough(const scalar KsPlus, const scalar Cs) const;

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("nutkRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutkRoughWallFunctionFvPatchScalarField
        //- onto a new patch
        nutkRoughWallFunctionFvPatchScalarField
        (
            const nutkRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutkRoughWallFunctionFvPatchScalarField
        (
            const nutkRoughWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutkRoughWallFunctionFvPatchScalarField
        (
            const nutkRoughWallFunctionFvPatchScalarField&,
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

        // Access

            //- Return Ks
            scalarField& Ks()
            {
                return Ks_;
            }

            //- Return Cs
            scalarField& Cs()
            {
                return Cs_;
            }


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
