/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::electrostaticDepositionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    The \c electrostaticDeposition is a boundary condition to
    calculate electric potential (\c V) on a given boundary
    based on film thickness (\c h) and film resistance (\c R) fields
    which are updated based on a given patch-normal current density
    field (\c jn), Coulombic efficiency and film resistivity.

    \f[
        j_n = - \sigma \nabla^\perp_p V = - \sigma (\vec{n}\cdot(\nabla V)_p)
    \f]

    \f[
        \frac{dh}{dt} = C_{eff} (j_n - j_{min})
    \f]

    \f[
        \frac{dR}{dt} = \rho \frac{dh}{dt} = \rho C_{eff} (j_n - j_{min})
    \f]

    \f[
        V_{film}^n = V_{film}^o + j_n R_\Delta
    \f]

    \f[
        V_{body} = j_n R_{body}
    \f]

    \f[
        V_p^n = V_i + V_{body} + V_{film}^n
    \f]

    where
    \vartable
      j_n     | Patch-normal current density        [A/m^2]
      V_p^n | Electric potential on film-fluid interface [volt = kg m^2/(A s^3)]
      V_p^o | Previous time-step electric potential on the interface [volt]
      V_{film} | Electric potential due to film resistance           [volt]
      V_{body} | Electric potential due to body resistance           [volt]
      V_i      | Initial electric potential                          [volt]
      R_\Delta| Film resistance (finite increment)  [ohm m^2 = kg m^4/(A^2 s^3)]
      R_{body} | Body resistance                    [ohm m^2 = kg m^4/(A^2 s^3)]
      \rho    | Isotropic film resistivity          [ohm m = kg m^3/(A^2 s^3)]
      h       | Film thickness                      [m]
      C_{eff} | Volumetric Coulombic efficiency     [m^3/(A s)]
      j_{min} | Minimum current density for deposition onset  [A/m^2]
      \sigma  | Isotropic conductivity of mixture   [S/m = A^2 s^3/(kg m^3)]
      \vec{n} | Patch-normal unit vector            [-]
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                    electrostaticDeposition;
        h                       <scalarField>;
        CoulombicEfficiency     <PatchFunction1>;
        resistivity             <PatchFunction1>;

        // Conditional mandatory entries

            // Option-1: single-phase
            sigma       <scalar>;

            // Option-2: multiphase
            phases
            {
                alpha.air
                {
                    sigma       <scalar>;
                }
                alpha.water
                {
                    sigma       <scalar>;
                }
                alpha.mercury
                {
                    sigma       <scalar>;
                }
                ...
            }

        // Optional entries
        jMin                    <scalar>;
        qMin                    <scalar>;
        Rbody                   <scalar>;
        Vi                      <scalar>;
        Vanode                  <scalar>;
        qCumulative             <scalarField>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type   | Reqd | Deflt
      type     | Type name: electrostaticDeposition  | word   | yes  | -
      h        | Film thickness                 | scalarField | yes  | -
      CoulombicEfficiency | Coulombic efficiency             <!--
             -->                   | PatchFunction1\<scalar\> | yes  | -
      resistivity         | Isotropic film resistivity       <!--
             -->                   | PatchFunction1\<scalar\> | yes  | -
      sigma  | Isotropic electrical conductivity of phase   | scalar | yes | -
      jMin   | Minimum current density for deposition onset | scalar | no  | 0
      qMin   | Minimum accumulative specific charge for deposition onset <!--
             -->                                             | scalar | no | 0
      Rbody  | Resistance due to main body and/or            <!--
             --> pretreatment layers                         | scalar | no | 0
      Vi     | Initial electric potential                    | scalar | no | 0
      Vanode | Anode electric potential                 | scalar | no | GREAT
      qCumulative  | Accumulative specific charge [A s/m^2]  <!--
             -->                                        | scalarField | no | 0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink

Note
  - Depletion or abrasion of material due to negative current is not allowed.
  - When accumulative specific charge (\c qCumulative) is less than minimum
  accumulative specific charge (\c qMin), no deposition occurs.
  - Boundary-condition updates are not allowed during outer corrections
  to prevent spurious accumulation of film thickness.
  - \c resistivity, \c jMin, \c qMin and \c Rbody are always non-negative.

SourceFiles
    electrostaticDepositionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef electrostaticDepositionFvPatchScalarField_H
#define electrostaticDepositionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Enum.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class electrostaticDepositionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class electrostaticDepositionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Film thickness
        scalarField h_;

        //- Accumulative specific charge
        scalarField qcum_;

        //- Electric potential due to film resistance
        scalarField Vfilm_;

        //- Coulombic efficiency
        autoPtr<PatchFunction1<scalar>> Ceffptr_;

        //- Isotropic film resistivity
        autoPtr<PatchFunction1<scalar>> rptr_;

        //- Minimum current density for the deposition onset
        scalar jMin_;

        //- Minimum accumulative specific charge for the deposition onset
        scalar qMin_;

        //- Resistance due to main body and/or pretreatment layers
        scalar Rbody_;

        //- Initial electric potential
        scalar Vi_;

        //- Anode electric potential
        scalar Vanode_;

        // Isotropic electrical conductivitiy properties

            //- Dictionary of phase data
            dictionary phasesDict_;

            //- List of phase names
            wordList phaseNames_;

            //- Unallocated list of phase fields
            UPtrList<volScalarField> phases_;

            //- List of isotropic electrical conductivity of phases
            PtrList<dimensionedScalar> sigmas_;

            //- Isotropic electrical conductivity of a single phase
            dimensionedScalar sigma_;

        //- Time index - used to prevent film accumulation during outer iters
        label timei_;

        //- Master patch ID
        mutable label master_;


    // Private Member Functions

        //- Return non-const access to an electrostaticDeposition patch
        electrostaticDepositionFvPatchScalarField& eVPatch
        (
            const label patchi
        ) const;

        //- Return non-const access to the master patch ID
        label& master() noexcept
        {
            return master_;
        }

        //- Set master patch ID
        void setMaster() const;

        //- Round scalars of a given scalar field to dcml points decimal
        void round(scalarField& fld, const scalar dcml=1e8) const;

        //- Write film thickness field to facilitate postprocessing
        void writeFilmFields() const;


public:

    //- Runtime type information
    TypeName("electrostaticDeposition");


    // Constructors

        //- Construct from patch and internal field
        electrostaticDepositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electrostaticDepositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- electrostaticDepositionFvPatchScalarField onto a new patch
        electrostaticDepositionFvPatchScalarField
        (
            const electrostaticDepositionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electrostaticDepositionFvPatchScalarField
        (
            const electrostaticDepositionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        electrostaticDepositionFvPatchScalarField
        (
            const electrostaticDepositionFvPatchScalarField&,
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

            //- Return const access to film thickness patch field
            const scalarField& h() const noexcept
            {
                return h_;
            }


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation

            //- Return the isotropic electrical conductivity field of mixture
            tmp<scalarField> sigma() const;

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
