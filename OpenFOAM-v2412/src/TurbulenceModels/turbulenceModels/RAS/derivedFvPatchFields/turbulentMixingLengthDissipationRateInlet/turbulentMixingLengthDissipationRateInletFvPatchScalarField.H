/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::turbulentMixingLengthDissipationRateInletFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides an inlet condition for turbulent kinetic
    energy dissipation rate, i.e. \c epsilon, based on a specified mixing
    length.  The patch values are calculated using:

        \f[
            \epsilon_p = \frac{C_{\mu}^{0.75} k^{1.5}}{L}
        \f]

    where
    \vartable
      \epsilon_p | Patch epsilon values     [m2/s3]
      C_\mu      | Empirical model constant retrived from turbulence model
      k          | Turbulent kinetic energy [m2/s2]
      L          | Mixing length scale      [m]
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            turbulentMixingLengthDissipationRateInlet;

        // Mandatory entries (runtime modifiable)
        mixingLength    0.005;

        // Optional entries (runtime modifiable)
        Cmu             0.09;
        k               k;
        phi             phi;

        // Placeholder
        value           uniform 200;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                  | Type   | Req'd | Dflt
      mixingLength | Mixing length scale [m]      | scalar |  yes  | -
      Cmu          | Empirical model constant     | scalar |  no   | 0.09
      phi          | Name of flux field           | word   |  no   | phi
      k            | Name of turbulent kinetic energy field | word | no | k
    \endtable

Note
    - The boundary condition is derived from \c inletOutlet condition.
      Therefore, in the event of reverse flow, a zero-gradient condition
      is applied.
    - The order of precedence to input the empirical model constant \c Cmu is:
      turbulence model, boundary condition dictionary, and default value=0.09.
    - The empirical model constant \c Cmu is not a spatiotemporal variant field.
      Therefore, the use of the boundary condition may not be fully consistent
      with the turbulence models where \c Cmu is a variant field, such as
      \c realizableKE closure model in this respect. Nevertheless, workflow
      observations suggest that the matter poses no importance.

See also
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentMixingLengthDissipationRateInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMixingLengthDissipationRateInletFvPatchScalarField_H
#define turbulentMixingLengthDissipationRateInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class turbulentMixingLengthDissipationRateInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingLengthDissipationRateInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Name of the turbulent kinetic energy field
        word kName_;

        //- Mixing length scale
        scalar mixingLength_;

        //- Empirical model constant
        scalar Cmu_;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthDissipationRateInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentMixingLengthDissipationRateInletFvPatchScalarField
        //- onto a new patch
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&,
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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
