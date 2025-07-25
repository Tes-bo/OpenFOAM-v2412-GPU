/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::atmTurbulentHeatFluxTemperatureFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a fixed heat constraint on temperature
    (i.e. \c T) to specify temperature gradient through an input heat source
    which can either be specified as absolute power [W], or as flux [W/m2].

    Required fields:
    \verbatim
      nut    | Turbulent viscosity                   [m2/s]
      k      | Turbulent kinetic energy              [m2/s2]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            atmTurbulentHeatFluxTemperature;
        heatSource      flux;
        alphaEff        alphaEff;

        // Mandatory entries (runtime modifiable)
        Cp0             0.001;
        q               uniform 10;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                         | Type | Reqd  | Dflt
      heatSource   | Heat source type                    | word |  yes  | -
      alphaEff     | Name of turbulent thermal diff. field [kg/m/s] <!--
                                                     --> | word | yes   | -
      Cp0          | Specific heat capacity [m2/s2/K] <!--
                                     --> | Function1<scalar>  | yes | -
      q            | Heat source value [W (power) or W/m2 (flux)] <!--
                                     --> | PatchFunction1<scalar> | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link fixedGradientFvPatchScalarField.H \endlink
     - \link Function1.H \endlink
     - \link PatchFunction1.H \endlink

    Options for the \c heatSource entry:
    \verbatim
      power  | Absolute power heat source [W]
      flux   | Flux heat source           [W/m2]
    \endverbatim

SourceFiles
    atmTurbulentHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmTurbulentHeatFluxTemperatureFvPatchScalarFields_H
#define atmTurbulentHeatFluxTemperatureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Function1.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmTurbulentHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmTurbulentHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private Enumerations

        //- Options for the heat source input modes
        enum heatSourceType
        {
            POWER,    //!< "Absolute power heat source"
            FLUX      //!< "Flux heat source"
        };

        //- Names for heatSourceType
        static const Enum<heatSourceType> heatSourceTypeNames;


    // Private Data

        //- Heat source type
        const enum heatSourceType heatSource_;

        //- Name of effective thermal diffusivity field [kg/m/s]
        word alphaEffName_;

        //- Specific heat capacity [m2/s2/K]
        autoPtr<Function1<scalar>> Cp0_;

        //- Heat power [W] or flux [W/m2]
        //  Divided by density, rho, if used in kinematic form
        autoPtr<PatchFunction1<scalar>> q_;

public:

    //- Runtime type information
    TypeName("atmTurbulentHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        atmTurbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmTurbulentHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmTurbulentHeatFluxTemperatureFvPatchScalarField onto
        //- a new patch
        atmTurbulentHeatFluxTemperatureFvPatchScalarField
        (
            const atmTurbulentHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmTurbulentHeatFluxTemperatureFvPatchScalarField
        (
            const atmTurbulentHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmTurbulentHeatFluxTemperatureFvPatchScalarField
        (
            const atmTurbulentHeatFluxTemperatureFvPatchScalarField&,
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


        // Evaluation

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
