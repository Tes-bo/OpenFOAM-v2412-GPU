/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

\defgroup grpLagrangianIntermediate Lagrangian particle models
@{
    This group contains Lagrangian modelling available in the 'intermediate' library
@}

\defgroup grpLagrangianIntermediateClouds Clouds
@{
    \ingroup grpLagrangianIntermediate
    This group contains Lagrangian clouds
@}

\defgroup grpLagrangianIntermediateParcels Parcels
@{
    \ingroup grpLagrangianIntermediate
    This group contains Lagrangian parcels
@}

// Submodels
\defgroup grpLagrangianIntermediateSubModels Submodels
@{
    \ingroup grpLagrangianIntermediate
    This group contains Lagrangian parcel submodels
@}

\defgroup grpLagrangianIntermediateKinematicSubModels Kinematic
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian kinematic parcel submodels
@}

\defgroup grpLagrangianIntermediateThermoSubModels Thermodynamic
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian thermo parcel submodels
@}

\defgroup grpLagrangianIntermediateReactingSubModels Reacting
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian reacting parcel submodels
@}

\defgroup grpLagrangianIntermediateReactingMultiphaseSubModels Reacting multiphase
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian reacting multiphase parcel submodels
@}

\defgroup grpLagrangianIntermediateMPPICSubModels MP-PIC
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian MP-PIC parcel submodels
@}


// Kinematic parcel sub models
\defgroup grpLagrangianIntermediateCollisionSubModels Collision
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle collision submodes
@}

\defgroup grpLagrangianIntermediateDispersionSubModels Dispersion
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle dispersion submodes
@}

\defgroup grpLagrangianIntermediateInjectionSubModels Injection
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle injection submodes
@}

\defgroup grpLagrangianIntermediateForceSubModels Forces
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle force submodels
@}

\defgroup grpLagrangianIntermediatePatchInteractionSubModels Patch interaction
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle patch interaction submodels
@}

\defgroup grpLagrangianIntermediateSurfaceFilmSubModels Surface film
@{
    \ingroup grpLagrangianIntermediateKinematicSubModels
    This group contains Lagrangian particle patch interaction submodels
@}


// Thermo parcel submodels
\defgroup grpLagrangianIntermediateHeatTransferSubModels Heat transfer
@{
    \ingroup grpLagrangianIntermediateThermoSubModels
    This group contains Lagrangian particle heat transfer submodels
@}


// Reacting parcel submodels
\defgroup grpLagrangianIntermediateCompositionSubModels Composition
@{
    \ingroup grpLagrangianIntermediateReactingSubModels
    This group contains Lagrangian particle composition submodels
@}

\defgroup grpLagrangianIntermediatePhaseChangeSubModels Phase change
@{
    \ingroup grpLagrangianIntermediateReactingSubModels
    This group contains Lagrangian particle phase change submodels
@}

// Spray parcel submodels (Spray derived from Reacting)
\defgroup grpLagrangianIntermediateAtomizationSubModels Atomization
@{
    \ingroup grpLagrangianIntermediateReactingSubModels
    This group contains Lagrangian particle atomization submodels
@}

\defgroup grpLagrangianIntermediateBreakupSubModels Breakup
@{
    \ingroup grpLagrangianIntermediateReactingSubModels
    This group contains Lagrangian particle breakup submodels
@}


// Reacting multiphase parcel submodels
\defgroup grpLagrangianIntermediateDevolatilisationSubModels Devolatilisation
@{
    \ingroup grpLagrangianIntermediateReactingMultiphaseSubModels
    This group contains Lagrangian particle devolatilisation submodels
@}

\defgroup grpLagrangianIntermediateSurfaceReactionSubModels Surface reaction
@{
    \ingroup grpLagrangianIntermediateReactingMultiphaseSubModels
    This group contains Lagrangian particle surface reaction submodels
@}


// MP-PIC parcel submodels
\defgroup grpLagrangianIntermediateMPPICAveragingMethods Averaging
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle averaging methods
@}

\defgroup grpLagrangianIntermediateMPPICCorrectionLimitingMethods Correction limiting
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle correction limiting methods
@}

\defgroup grpLagrangianIntermediateMPPICDampingSubModels Damping
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle damping submodels
@}

\defgroup grpLagrangianIntermediateMPPICIsotropySubModels Isotropy
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle isotropy submodels
@}

\defgroup grpLagrangianIntermediateMPPICPackingSubModels Packing
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle packing submodels
@}

\defgroup grpLagrangianIntermediateMPPICParticleStressSubModels Particle stress
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle stress submodels
@}

\defgroup grpLagrangianIntermediateMPPICTimeScaleSubModels Time scale
@{
    \ingroup grpLagrangianIntermediateMPPICSubModels
    This group contains Lagrangian MP-PIC particle time scale submodels
@}


// Function objects
\defgroup grpLagrangianIntermediateFunctionObjects Function objects
@{
    \ingroup grpLagrangianIntermediateSubModels
    This group contains Lagrangian function objects
@}


\*---------------------------------------------------------------------------*/
