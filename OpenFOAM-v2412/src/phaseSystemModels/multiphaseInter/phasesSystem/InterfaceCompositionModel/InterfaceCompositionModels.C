/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"
#include "Boussinesq.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "solidThermo.H"
#include "heSolidThermo.H"
#include "solidThermoPhysicsTypes.H"

#include "kineticGasEvaporation.H"
#include "Lee.H"
#include "interfaceHeatResistance.H"
#include "interfaceOxideRate.H"
#include "diffusionGasEvaporation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > constRhoHThermoPhysics;


    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    Boussinesq<specie>
                >,
                sensibleEnthalpy
            >
        > BoussinesqFluidEThermoPhysics;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace meltingEvaporationModels;

    //NOTE: First thermo (from) and second otherThermo (to)

    // kineticGasEvaporation model definitions

        // From pure liquid (rhoConst) to a multi-component gas incomp phase
        makeInterfaceContSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );

        // From pure liquid (BoussinesqFluid) to a multi-component gas incomp
        // phase
        makeInterfaceContSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );


        // From pure liquid (rhoConst) to pure gas (incompressible ideal gas)
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics
        );

        // From pure liquid (const rho) to pure gas (rhoConst) gas
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );


        // From pure liquid (Boussinesq) to pure gas (incompressible ideal gas)
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics
        );

        // From pure liquid (Boussinesq) to pure gas (rho const)
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );


    // Lee model definitions

        // From pure phase (rho const) to phase (rho const)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure phase (rho const) to phase (Boussinesq)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );


        // From pure solid phase (const) to phase (Boussinesq)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );

        // From pure solid phase (const) to phase (rho const)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure solid phase (const) to phase (tabulated)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics
        );

        // From pure solid phase (const) to phase (poly)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics
        );

        // From pure solid phase (poly) to flow phase (poly)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPolyTranspPolyIcoSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics
        );

        // From pure solid phase (poly) to flow phase (tabulated)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPolyTranspPolyIcoSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics
        );

        // From pure solid tabulated to pure tabulated flow
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedPolyIcoSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics
        );

        // From pure phase (exp-Transp, hPower solidThermo) to phase (ico-rho)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPowerSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics
        );


        // From pure phase (const rho) to multi phase (incomp ideal gas)
        makeInterfaceContSpecieMixtureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );


        // From pure phase (Boussinesq) to phase (solidThermo)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        // From pure phase (rho const) to phase (solidThermo)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        //From pure phase (poly) to solid phase (exp)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPowerSolidThermoPhysics
        );

        //From pure phase (poly) to solid phase (poly)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPolyTranspPolyIcoSolidThermoPhysics
        );

        //From pure phase (poly) to solid phase (const)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        //From pure phase (tabulated) to solid phase (poly)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPolyTranspPolyIcoSolidThermoPhysics
        );

        //From pure phase (tabulated) to solid phase (const)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedPolyIcoSolidThermoPhysics
        );

        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedPolyIcoSolidThermoPhysics
        );

        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedPolyIcoSolidThermoPhysics
        );

        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics
        );

        // From pure phase (tabulated) to solid phase (icoTabulated)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedIcoTabulatedSolidThermoPhysics
        );

        // From solid phase (icoTabulated) to pure phase (tabulated)
        makeInterfacePureType
        (
            Lee,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hTabulatedIcoTabulatedSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics
        );


        // interfaceHeatResistance model definitions

        // From pure phase (rho const) to phase (rho const)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure phase (rho const) to phase (Boussinesq)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );


        // From pure phase (solidThermo) to phase (Boussinesq)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );

        // From pure phase (solidThermo) to phase (rho const)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure phase (all-poly solidThermo) to phase (ico-rho)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPolyTranspPolyIcoSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics
        );

        // From pure phase (exp-Transp, hPower solidThermo) to phase (ico-rho)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPowerSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics
        );


        // From pure phase (const rho) to multi phase (incomp ideal gas)
        makeInterfaceContSpecieMixtureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );


        // From pure phase (Boussinesq) to phase (solidThermo)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        // From pure phase (rho const) to phase (solidThermo)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        //From pure liquid phase (ico-rho) to pure phase (exp-Transp, hPower solidThermo)
        makeInterfacePureType
        (
            interfaceHeatResistance,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hPowerSolidThermoPhysics
        );


        // interfaceOxideRate model definitions

        // From pure phase (tabulated) to solid phase (const)
        makeInterfacePureType
        (
            interfaceOxideRate,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        // From pure phase (rho const) to phase (rho const)
        makeInterfacePureType
        (
            interfaceOxideRate,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure phase (ico) to solid phase (const)
        makeInterfacePureType
        (
            interfaceOxideRate,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        // From pure phase (tabulated) to phase (rho const)
        makeInterfacePureType
        (
            interfaceOxideRate,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            tabulatedThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );

        // From pure phase (ico) to phase (rho const)
        makeInterfacePureType
        (
            interfaceOxideRate,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            icoPoly8HThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics
        );


    // diffusionGasEvaporation model definitions

        // From pure liquid (rhoConst) to a multi-component gas incomp phase
        makeInterfaceContSpecieMixtureType
        (
            diffusionGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );


        // From pure liquid (BoussinesqFluid) to a multi-component gas incomp
        // phase
        makeInterfaceContSpecieMixtureType
        (
            diffusionGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );
}


// ************************************************************************* //
