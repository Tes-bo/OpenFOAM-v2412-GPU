/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2022 OpenCFD Ltd.
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

#include "makeSolidThermo.H"
#include "solidThermo.H"
#include "pureZoneMixture.H"
#include "heSolidThermo.H"

#include "specie.H"
#include "rhoConst.H"
#include "icoPolynomial.H"
#include "icoTabulated.H"
#include "hConstThermo.H"
#include "hPowerThermo.H"
#include "hPolynomialThermo.H"
#include "hTabulatedThermo.H"
#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
#include "exponentialSolidTransport.H"
#include "polynomialSolidTransport.H"
#include "tabulatedSolidTransport.H"
#include "tabulatedAnIsoSolidTransport.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

// Note: pureZoneMixture copies mixture model for every evaluation of cell
//       so can become expensive for complex models (e.g. with tables)

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    constAnIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    exponentialSolidTransport,
    sensibleEnthalpy,
    hPowerThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    polynomialSolidTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    tabulatedSolidTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    icoPolynomial,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureZoneMixture,
    tabulatedSolidTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    icoTabulated,
    specie
);

//- TBD. Needs clone
//makeSolidThermo
//(
//    solidThermo,
//    heSolidThermo,
//    pureZoneMixture,
//    tabulatedAnIsoSolidTransport,
//    sensibleEnthalpy,
//    hTabulatedThermo,
//    icoPolynomial,
//    specie
//);

//makeSolidThermoPhysicsType
//(
//    solidThermo,
//    heSolidThermo,
//    pureZoneMixture,
//    hTransportThermoPoly8SolidThermoPhysics
//);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
