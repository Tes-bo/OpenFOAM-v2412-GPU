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
    Foam::nutURoughWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function on the turbulent
    viscosity (i.e. \c nut) based on velocity (i.e. \c U) for low- and
    high-Reynolds number applications for rough walls.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                nutURoughWallFunction;
        roughnessHeight     1e-5;
        roughnessConstant   0.5;
        roughnessFactor     1;

        // Optional entries
        maxIter             10;
        tolerance           0.0001;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd   | Deflt
      type      | Type name: nutURoughWallFunction    | word | yes    | -
      roughnessHeight   | Roughness height            | scalar | yes  | -
      roughnessConstant | Roughness constant          | scalar | yes  | -
      roughnessFactor   | Scaling factor              | scalar | yes  | -
      maxIter   | Number of Newton-Raphson iterations | label  | no   | 10
      tolerance | Convergence tolerance               | scalar | no   | 0.0001
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink

Note
    - Suffers from non-exact restart since \c correctNut() (called through
    \c turbulence->validate) returns a slightly different value every time
    it is called.
    See \link nutUSpaldingWallFunctionFvPatchScalarField.C \endlink.

SourceFiles
    nutURoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutURoughWallFunctionFvPatchScalarField_H
#define nutURoughWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutURoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutURoughWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private Data

        // Roughness model parameters

            //- Height
            scalar roughnessHeight_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFactor_;


        //- Max iterations in calcNut
        const label maxIter_;

        //- Convergence tolerance
        const scalar tolerance_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate yPlus
        tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutURoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutURoughWallFunctionFvPatchScalarField
        //- onto a new patch
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
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

            //- Return the roughness height
            scalar roughnessHeight() const noexcept
            {
                return roughnessHeight_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight() noexcept
            {
                return roughnessHeight_;
            }

            //- Return the roughness constant scale
            scalar roughnessConstant() const noexcept
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant() noexcept
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            scalar roughnessFactor() const noexcept
            {
                return roughnessFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //- adjustment
            scalar& roughnessFactor() noexcept
            {
                return roughnessFactor_;
            }


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
