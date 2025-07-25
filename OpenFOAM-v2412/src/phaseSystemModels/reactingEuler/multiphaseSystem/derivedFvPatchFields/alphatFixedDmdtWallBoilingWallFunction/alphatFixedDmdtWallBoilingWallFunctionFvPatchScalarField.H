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
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField

Description
    A simple \c alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField with
    a fixed volumetric phase-change mass flux.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patch>
    {
        // Mandatory entries
        type         compressible::alphatFixedDmdtWallBoilingWallFunction;
        vaporPhase   <word>;

        // Optional entries
        relax        <scalar>;
        fixedDmdt    <scalar>;
        L            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      type      | Type name:                             <!--
      -->  compressible::alphatFixedDmdtWallBoilingWallFunction <!--
      -->                                                 | word | yes  | -
      vaporPhase | Name of the vapor phase                | word | yes  | -
      relax      | Relaxation factor for dmdt           | scalar | no   | 1.0
      fixedDmdt  | Volumetric phase-change mass flux in near wall cells <!--
      -->                                               | scalar | no   | 0.0
      L          | Latent heat                          | scalar | no   | 0.0
    \endtable

    The inherited entries are elaborated in:
     -\link alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H\endlink

See also
    Foam::compressible::
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField

SourceFiles
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H
#define alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
{
    // Private Data

        //- Name of the vapor phase
        word vaporPhaseName_;

        //- Relaxation factor for dmdt
        scalar relax_;

        //- Volumetric phase-change mass flux in near wall cells
        scalar fixedDmdt_;

        //- Latent heat
        scalar L_;


public:

    //- Runtime type information
    TypeName("compressible::alphatFixedDmdtWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        //- onto a new patch
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&,
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

        //- Is there phase change mass transfer for this phasePair
        virtual bool activePhasePair(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& dmdt(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& mDotL(const phasePairKey&) const;


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
