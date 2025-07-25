/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

Typedefs
    Foam::thermoPhysicsTypes

Description
    Type definitions for thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef thermoPhysicsTypes_H
#define thermoPhysicsTypes_H

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "incompressiblePerfectGas.H"
#include "rPolynomial.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "tabulatedTransport.H"
#include "hTabulatedThermo.H"
#include "icoTabulated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on sensibleEnthalpy
    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constGasHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > gasHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                PengRobinsonGas<specie>
            >,
            sensibleEnthalpy
        >
    > PengRobinsonGasHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleGasHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > incompressibleGasHThermoPhysics;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8HThermoPhysics;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >,
        8
    > icoPoly8TranspJanafHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constFluidHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                rPolynomial<specie>
            >,
            sensibleEnthalpy
        >
    >
    constrPolFluidHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    >
    constAdiabaticFluidHThermoPhysics;

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
    >
    constHThermoPhysics;

    typedef
    tabulatedTransport
    <
        species::thermo
        <
            hTabulatedThermo
            <
                icoTabulated<specie>
            >,
            sensibleEnthalpy
        >
    >
    tabulatedThermoPhysics;


    // thermo physics types based on sensibleInternalEnergy
    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constGasEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > gasEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                PengRobinsonGas<specie>
            >,
            sensibleInternalEnergy
        >
    > PengRobinsonGasEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleGasEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                incompressiblePerfectGas<specie>
            >,
            sensibleInternalEnergy
        >
    > incompressibleGasEThermoPhysics;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            hPolynomialThermo
            <
                icoPolynomial<specie, 8>,
                8
            >,
            sensibleInternalEnergy
        >,
        8
    > icoPoly8EThermoPhysics;

    typedef
    polynomialTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleInternalEnergy
        >,
        8
    > icoPoly8TranspJanafEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constrPolFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constAdiabaticFluidEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            eConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    >
    constEThermoPhysics;

    typedef
    tabulatedTransport
    <
        species::thermo
        <
            hTabulatedThermo
            <
                icoTabulated<specie>
            >,
            sensibleInternalEnergy
        >
    >
    tabulatedEThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
