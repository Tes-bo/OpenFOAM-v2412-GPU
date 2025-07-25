/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
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
    Foam::atmNutkWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    viscosity (i.e. \c nut) based on the turbulent kinetic energy (i.e. \c k)
    for atmospheric boundary layer modelling.  It is designed to be used in
    conjunction with the \c atmBoundaryLayerInletVelocity boundary condition.

    The governing equation of the boundary condition:

        \f[
            u = \frac{u^*}{\kappa} ln \left(\frac{z + z_0}{z_0}\right)
        \f]

    where
    \vartable
      u^*      | Friction velocity
      \kappa   | von Kármán constant
      z_0      | Surface roughness length       [m]
      z        | Ground-normal coordinate
    \endvartable

    Required fields:
    \verbatim
      nut    | Turbulent viscosity              [m2/s]
      k      | Turbulent kinetic energy         [m2/s2]
    \endverbatim

    References:
    \verbatim
        Theoretical expressions (tag:HW):
            Hargreaves, D. M., & Wright, N. G. (2007).
            On the use of the k–ε model in commercial CFD software to model
            the neutral atmospheric boundary layer.
            J. of wind engineering and industrial aerodynamics, 95(5), 355-369.
            DOI:10.1016/j.jweia.2006.08.002
    \endverbatim

    Required fields:
    \verbatim
      nut    | Turbulent viscosity              [m2/s]
      k      | Turbulent kinetic energy         [m2/s2]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmNutkWallFunction;
        z0              <PatchFunction1<scalar>>;

        // Optional entries
        boundNut        <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd    | Deflt
      type     | Type name: atmNutkWallFunction   | word   | yes     | -
      z0       | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
      boundNut | Flag to zero-bound nut near wall | bool   | no      | false
    \endtable

    The inherited entries are elaborated in:
     - \link nutkWallFunctionFvPatchScalarField.H \endlink
     - \link PatchFunction1.H \endlink

Note
    - \c boundNut entry is set \c false for backward compatibility reasons.
    - \c nutkAtmRoughWallFunction was renamed to \c atmNutkWallFunction.

SourceFiles
    atmNutkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmNutkWallFunctionFvPatchScalarField_H
#define atmNutkWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class atmNutkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmNutkWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
    // Private Data

        //- Flag to zero-bound nut to prevent negative nut
        //- at the wall arising from negative heat fluxes
        const bool boundNut_;

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
    TypeName("atmNutkWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmNutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmNutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmNutkWallFunctionFvPatchScalarField
        //- onto a new patch
        atmNutkWallFunctionFvPatchScalarField
        (
            const atmNutkWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmNutkWallFunctionFvPatchScalarField
        (
            const atmNutkWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmNutkWallFunctionFvPatchScalarField
        (
            const atmNutkWallFunctionFvPatchScalarField&,
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
