/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::solidThermoPhysicsTypes

Description
    Type definitions for solid-thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef solidThermoPhysicsTypes_H
#define solidThermoPhysicsTypes_H

#include "specie.H"
#include "rhoConst.H"
#include "icoPolynomial.H"
#include "icoTabulated.H"
#include "hConstThermo.H"
#include "hPolynomialThermo.H"
#include "hPowerThermo.H"
#include "hTabulatedThermo.H"
#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
#include "exponentialSolidTransport.H"
#include "polynomialSolidTransport.H"
#include "tabulatedSolidTransport.H"

#include "sensibleEnthalpy.H"
#include "thermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        constIsoSolidTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > hConstSolidThermoPhysics;

    typedef
        exponentialSolidTransport
        <
            species::thermo
            <
                hPowerThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > hPowerSolidThermoPhysics;

    typedef
        polynomialSolidTransport
        <
            species::thermo
            <
                hPolynomialThermo
                <
                    rhoConst<specie>,
                    8
                >,
                sensibleEnthalpy
            >,
            8
        > hTransportThermoPoly8SolidThermoPhysics;

     typedef
        constIsoSolidTransport
        <
            species::thermo
            <
                hPowerThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > hExpKappaConstSolidThermoPhysics;

    typedef
        polynomialSolidTransport
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
        > hPolyTranspPolyIcoSolidThermoPhysics;

     typedef
        tabulatedSolidTransport
        <
            species::thermo
            <
                hTabulatedThermo
                <
                    icoPolynomial<specie, 8>
                >,
                sensibleEnthalpy
            >
        > hTabulatedPolyIcoSolidThermoPhysics;

    typedef
        tabulatedSolidTransport
        <
            species::thermo
            <
                hTabulatedThermo
                <
                    icoTabulated<specie>
                >,
                sensibleEnthalpy
            >
        > hTabulatedIcoTabulatedSolidThermoPhysics;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
