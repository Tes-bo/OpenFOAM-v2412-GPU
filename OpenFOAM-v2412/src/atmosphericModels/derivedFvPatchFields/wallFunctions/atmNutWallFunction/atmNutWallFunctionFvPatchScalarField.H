/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 CENER
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
    Foam::atmNutWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    viscosity (i.e. \c nut) based on the turbulent kinetic energy (i.e. \c k)
    and velocity (i.e. \c U) for atmospheric boundary layer modelling.

    The governing equation of the boundary condition:

        \f[
            \tau_w = {U^*_u} {U^*_k}
        \f]

    with

        \f[
            {U^*_u} = \frac{\kappa U_w}{ln(z_p / z_0)}
        \f]

        \f[
            {U^*_k} = C_{\mu}^{1/4} \sqrt{k}
        \f]

    where
    \vartable
        \tau_w  | wall shear stress
        U^*_u   | local friction velocity based on near-ground velocity
        U^*_k   | local friction velocity based on near-ground k
        \kappa  | von Kármán constant
        U_w     | near-ground velocity
        z_p     | vertical coordinate
        z_0     | surface roughness length [m]
        C_mu    | empirical model constant
        k       | turbulent kinetic energy
    \endvartable

    References:
    \verbatim
        Theoretical expressions (tags:RH, SBJM, SM):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

            Sørensen, N. N., Bechmann, A., Johansen, J., Myllerup, L.,
            Botha, P., Vinther, S., & Nielsen, B. S. (2007).
            Identification of severe wind conditions using
            a Reynolds Averaged Navier-Stokes solver.
            In Journal of Physics: Conference
            series (Vol. 75, No. 1, p. 012053).
            DOI:10.1088/1742-6596/75/1/012053

            Sumner, J., & Masson, C. (2012).
            k−ε simulations of the neutral atmospheric boundary layer:
            analysis and correction of discretization errors on practical grids.
            International journal for numerical
            methods in fluids, 70(6), 724-741.
            DOI:10.1002/fld.2709
    \endverbatim

    Required fields:
    \verbatim
      nut     | Turbulent viscosity      [m2/s]
      k       | Turbulent kinetic energy [m2/s2]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmNutWallFunction;
        z0Min           <scalar>;
        z0              <PatchFunction1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                       | Type   | Reqd  | Deflt
      type     | Type name: nutAtmWallFunction     | word   | yes   | -
      z0Min | Minimum surface roughness length [m] | scalar | yes   | -
      z0    | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link nutkWallFunctionFvPatchScalarField.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    atmNutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmNutWallFunctionFvPatchScalarField_H
#define atmNutWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class atmNutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmNutWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{

    // Private Data

        //- Minimum surface roughness length [m]
        const scalar z0Min_;

        //- Surface roughness length field [m]
        autoPtr<PatchFunction1<scalar>> z0_;


protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("atmNutWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmNutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmNutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmNutWallFunctionFvPatchScalarField
        //- onto a new patch
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&,
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
