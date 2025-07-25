/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::atmOmegaWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the specific
    dissipation rate (i.e. \c omega) and the turbulent kinetic energy
    production contribution (i.e. \c G) for atmospheric boundary
    layer modelling.

    References:
    \verbatim
        Theoretical expressions (tags:PGVB, B):
            Parente, A., Gorlé, C., Van Beeck, J., & Benocci, C. (2011).
            Improved k–ε model and wall function formulation
            for the RANS simulation of ABL flows.
            J. of wind engineering and industrial aerodynamics, 99(4), 267-278.
            DOI:10.1016/j.jweia.2010.12.017

            Bredberg, J. (2000).
            On the wall boundary condition for turbulence models.
            Chalmers University of Technology, Depart. of Thermo and Fluid Dyn.
            Internal Report 00/4. Sweden: Göteborg.
    \endverbatim

    Required fields:
    \verbatim
      omega    | Specific dissipation rate    [1/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmOmegaWallFunction;
        z0              <PatchFunction1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd   | Deflt
      type     | Type name: atmOmegaWallFunction  | word   | yes    | -
      z0    | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link omegaWallFunctionFvPatchScalarField.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    atmOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmOmegaWallFunctionFvPatchScalarField_H
#define atmOmegaWallFunctionFvPatchScalarField_H

#include "omegaWallFunctionFvPatchScalarField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class atmOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmOmegaWallFunctionFvPatchScalarField
:
    public omegaWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Surface roughness length field [m]
        autoPtr<PatchFunction1<scalar>> z0_;


    // Protected Member Functions

        //- Calculate the omega and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("atmOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmOmegaWallFunctionFvPatchScalarField
        //- onto a new patch
        atmOmegaWallFunctionFvPatchScalarField
        (
            const atmOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmOmegaWallFunctionFvPatchScalarField
        (
            const atmOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmOmegaWallFunctionFvPatchScalarField
        (
            const atmOmegaWallFunctionFvPatchScalarField&,
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
    virtual ~atmOmegaWallFunctionFvPatchScalarField() = default;


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
