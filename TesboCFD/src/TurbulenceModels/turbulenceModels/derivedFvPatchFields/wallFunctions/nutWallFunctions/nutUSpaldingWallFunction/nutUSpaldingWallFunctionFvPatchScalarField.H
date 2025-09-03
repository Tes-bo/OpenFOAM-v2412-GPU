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
    Foam::nutUSpaldingWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the turbulent
    viscosity (i.e. \c nut) based on velocity (i.e. \c U). Using Spalding's
    law gives a continuous \c nut profile to the wall.

    \f[
        y^+ = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
            - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
    \f]

    where
    \vartable
      y^+     | Wall-normal height of a cell centre in wall units
      u^+     | Velocity at \f$y^+\f$ in wall units
      \kappa  | von Karman constant
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutUSpaldingWallFunction;

        // Optional entries
        maxIter         10;
        tolerance       0.0001;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type   | Reqd | Deflt
      type      | Type name: nutUSpaldingWallFunction | word   | yes  | -
      maxIter   | Number of Newton-Raphson iterations | label  | no   | 10
      tolerance | Convergence tolerance               | scalar | no   | 0.0001
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink

Note
  - Suffers from non-exact restart since \c correctNut() (called through
    \c turbulence->validate) returns a slightly different value every time
    it is called. This is since the seed for the Newton-Raphson iteration
    uses the current value of \c *this (\c =nut ).
  - This can be avoided by overriding the tolerance. This also switches on
    a pre-detection whether the current nut already satisfies the turbulence
    conditions and if so does not change it at all. This means that the nut
    only changes if it 'has really changed'. This probably should be used with
    a tight tolerance, to make sure to kick every iteration, e.g.
        maxIter     100;
        tolerance   1e-7;

SourceFiles
    nutUSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUSpaldingWallFunctionFvPatchScalarField_H
#define nutUSpaldingWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUSpaldingWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Max iterations in calcNut
        const label maxIter_;

        //- Convergence tolerance
        const scalar tolerance_;

        //- Uncomment in case of intrumentation
        //mutable uint64_t invocations_;
        //mutable uint64_t nontrivial_;
        //mutable uint64_t nonconvergence_;
        //mutable uint64_t iterations_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Calculate the friction velocity and number of iterations for
        //- convergence
        tmp<scalarField> calcUTau
        (
            const scalarField& magGradU,
            const label maxIter,
            scalarField& err
        ) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutUSpaldingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutUSpaldingWallFunctionFvPatchScalarField
        //- onto a new patch
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&,
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


    //- Destructor
    virtual ~nutUSpaldingWallFunctionFvPatchScalarField();


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
