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
    Foam::atmEpsilonWallFunctionFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent kinetic
    energy dissipation rate (i.e. \c epsilon) and the turbulent kinetic energy
    production contribution (i.e. \c G) for atmospheric boundary layer
    modelling.

    References:
    \verbatim
        Theoretical expressions (tags:PGVB, RH):
            Parente, A., Gorlé, C., Van Beeck, J., & Benocci, C. (2011).
            Improved k–ε model and wall function formulation
            for the RANS simulation of ABL flows.
            J. of wind engineering and industrial aerodynamics, 99(4), 267-278.
            DOI:10.1016/j.jweia.2010.12.017

            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8
    \endverbatim

    Required fields:
    \verbatim
      epsilon    | Turbulent kinetic energy dissipation rate    [m2/s3]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            atmEpsilonWallFunction;
        z0              <PatchFunction1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                       | Type   | Reqd     | Deflt
      type     | Type name: atmEpsilonWallFunction | word   | yes      | -
      z0       | Surface roughness length [m] | PatchFunction1<scalar> | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link epsilonWallFunctionFvPatchScalarField.H \endlink
      - \link PatchFunction1 \endlink

SourceFiles
    atmEpsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmEpsilonWallFunctionFvPatchScalarField_H
#define atmEpsilonWallFunctionFvPatchScalarField_H

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class atmEpsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmEpsilonWallFunctionFvPatchScalarField
:
    public epsilonWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Surface roughness length field [m]
        autoPtr<PatchFunction1<scalar>> z0_;


    // Protected Member Functions

        //- Calculate the epsilon and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& epsilon
        );

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("atmEpsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmEpsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmEpsilonWallFunctionFvPatchScalarField
        //- onto a new patch
        atmEpsilonWallFunctionFvPatchScalarField
        (
            const atmEpsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmEpsilonWallFunctionFvPatchScalarField
        (
            const atmEpsilonWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmEpsilonWallFunctionFvPatchScalarField
        (
            const atmEpsilonWallFunctionFvPatchScalarField&,
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
    virtual ~atmEpsilonWallFunctionFvPatchScalarField() = default;


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
