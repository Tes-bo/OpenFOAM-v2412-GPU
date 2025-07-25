/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::
    compressible::
    humidityTemperatureCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature to be used at the coupling
    interface between fluid solid regions.

    This boundary condition can operate in four modes:
    - \c constantMass: thermal inertia only
      - requires \c rho, \c thickness and \c cp
    - \c condensation: condensation only
      - when the wall temperature (Tw) is below the dew temperature (Tdew)
        condesation takes place and the resulting condensed mass is stored
        on the wall
    - \c evaporation: evaporation only
      - initial mass is vaporized when Tw is above the input vaporization
        temperature (Tvap).
    - \c condensationAndEvaporation : condensation and evaporation take place
      simultaneously.

    There is no mass flow on the wall, i.e. the mass condensed on a face
    remains on that face. It uses a 'lumped mass' model to include thermal
    inertia effects.

    It assumes a drop-wise type of condensation, whereby its heat transfer
    Nusselt number is calculated using:
    \f{eqnarray*}{
        51104 + 2044 (T - 273.15)   & T > 295 & T < 373 \\
        255510                      & T > 373 &
    \f}

    References:
    \verbatim
        Standard models (tag:BLID):
            Bergman, T. L., Lavine, A. S.,
            Incropera, F. P., & Dewitt, D. P. (2011).
            Fundamentals of heat and mass transfer.
            John Wiley & Sons. 7th Edition. Chapter 10.
            ISBN:9780470501979
    \endverbatim

    The mass transfer correlation used is:

    \f[ h_m = D_{ab} \frac{Sh}{L} \f]

    where:
    \vartable
        D_{ab} | mass vapour difussivity
        L      | characteristic length
        Sh     | Sherwood number
    \endvartable

    The Sherwood number is calculated using:

    \f{eqnarray*}{
            0.664 Re^\frac{1}{2} Sc^\frac{1}{3} & Re < 5.0E+05 \\
            0.037 Re^\frac{4}{5} Sc^\frac{1}{3} & Re > 5.0E+05
    \f}
    where:
    \vartable
        Re     | Reynolds number
        Sc     | Schmidt number
    \endvartable

    NOTE:
    - The correlation used to calculate Tdew is for water vapour.
    - A scalar transport equation for the carrier specie is required, e.g.
      supplied via a function object or in the main solver. This specie
      transports the vapour phase in the main ragion.
    - The boundary condition of this specie on the coupled wall must be
      fixedGradient in order to allow condensation or evaporation of the
      vapour in or out of this wall
    - Addition of extra layers in possible using thicknessLayers and
      kappaLayers


    Example usage:

    On the fluid side
    \verbatim
    myInterfacePatchName
    {
        type            thermalHumidityCoupledMixed;
        kappaMethod     fluidThermo;
        kappa           none;

        // Modes of operation: inert, condensation, vaporization, condEvap
        mode            condEvap;

        // Carrier species name
        specieName      H2O;

        // Carrier molecular weight
        carrierMolWeight           28.9;

        // Characteristic length of the wall
        L               0.1;

        // Vaporisation temperature
        Tvap            273;

        // Liquid properties for the condensed mass
        liquid
        {
            H2O
            {
                defaultCoeffs       yes;
            }
        }

        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);

        // thickness, density and cp required for inert and condensation
        // modes

        //thickness       uniform 0;
        //cp              uniform 0;
        //rho             uniform 0;

        value           $internalField;
    }
    \endverbatim

    On the solid side:
    \verbatim
    myInterfacePatchName
    {
        type            thermalInertiaMassTransferCoupledMixed;
        kappaMethod     solidThermo;
        kappa           none;
        value           uniform 260;
    }
    \endverbatim


SourceFiles
    humidityTemperatureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef humidityTemperatureCoupledMixedFvPatchScalarField_H
#define humidityTemperatureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "liquidProperties.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class humidityTemperatureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class humidityTemperatureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    //  Public enumeration

        //- Modes of mass transfer
        enum massTransferMode
        {
            mtConstantMass,
            mtCondensation,
            mtEvaporation,
            mtCondensationAndEvaporation
        };


private:

    // Private data

        static const Enum<massTransferMode> massModeTypeNames_;

        //- Operating mode
        massTransferMode mode_;


        // Field names

            //- Name of the pressure field
            const word pName_;

            //- Name of the velocity field
            const word UName_;

            //- Name of the density field
            const word rhoName_;

            //- Name of the dynamic viscosity field
            const word muName_;

            //- Name of temperature field on the neighbour region
            const word TnbrName_;

            //- Name of the radiative heat flux in the neighbour region
            const word qrNbrName_;

            //- Name of the radiative heat flux field
            const word qrName_;

            //- Name of the species on which the mass transferred (default H2O)
            const word specieName_;


        //- Liquid properties
        autoPtr<liquidProperties> liquid_;

        //- Liquid dictionary
        dictionary liquidDict_;

        //- Mass accumulated on faces
        scalarField mass_;

        //- Vaporization temperature
        scalar Tvap_;

        //- Cache myDelta
        scalarField myKDelta_;

        //- Phase change energy
        scalarField dmHfg_;

        //- Thermal inertia
        scalarField mpCpTp_;

        //- Average molecular weight for the carrier mixture in the gas phase
        scalar Mcomp_;

        //- Characteristic length scale
        scalar L_;

        //- Fluid side
        bool fluid_;

        //- Cp field for inert mode
        scalarField cp_;

        //- Thickness field for inert mode
        scalarField thickness_;

        //- Density field for inert mode
        scalarField rho_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;


    // Private members

        //- Calculation of Sh
        scalar Sh(const scalar Re, const scalar Sc) const;

        //- Calculation of htc from the condensed surface
        scalar htcCondensation(const scalar TSat, const scalar Re) const;

        //- Lookup (or create) thickness field for output
        volScalarField& thicknessField(const word&, const fvMesh&);


public:

    //- Runtime type information
    TypeName("humidityTemperatureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const
            humidityTemperatureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const humidityTemperatureCoupledMixedFvPatchScalarField&,
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


    // Member functions

            // Mapping functions

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


        //- Return myKDelta
        const scalarField myKDelta() const
        {
            return myKDelta_;
        }

        //- Return mpCpTp
        const scalarField mpCpTp() const
        {
            return mpCpTp_;
        }

        //- Return dmHfg
        const scalarField dmHfg() const
        {
            return dmHfg_;
        }

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
