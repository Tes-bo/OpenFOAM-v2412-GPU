/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "PhaseCompressibleTurbulenceModel.H"
#include "phaseModel.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Base models defined in compressibleTwoPhaseSystem
//
// Additional models only

// Typedefs
defineTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (phaseModelPhaseCompressibleTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, LES, Type)

// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

// #include "Stokes.H"


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSato.H"
makeRASModel(kOmegaSSTSato);

#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);

#include "LaheyKEpsilon.H"
makeRASModel(LaheyKEpsilon);

#include "continuousGasKEpsilon.H"
makeRASModel(continuousGasKEpsilon);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "SmagorinskyZhang.H"
makeLESModel(SmagorinskyZhang);

#include "NicenoKEqn.H"
makeLESModel(NicenoKEqn);

#include "continuousGasKEqn.H"
makeLESModel(continuousGasKEqn);


// -------------------------------------------------------------------------- //
// Additional models
// -------------------------------------------------------------------------- //

#include "kineticTheoryModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, kineticTheoryModel);

#include "phasePressureModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, phasePressureModel);


// ************************************************************************* //
