<!------------------------------------------------------------------------- -->

# Overview

"By setting appropriate profiles for wind velocity and the turbulence
quantities at the inlet, it is often assumed that the boundary layer will
be maintained up to the buildings or obstructions in the flow." (HW:p. 355).
However, it was quantified by (HW:p. 355) that "even in the absence of
obstructions, ..., the velocity and turbulence profiles decay along the
fetch" (HW:p. 355). It was shown by (HW:p. 355) that a set of modifications
were required to maintain a neutral atmospheric boundary layer throughout
an empty and long computational domain of a RANS computation.

## Aim

- Verification of the atmospheric boundary-layer boundary conditions in terms
of the maintenance of inlet quantities downstream within a RANS computation:
  - atmBoundaryLayerInletVelocity
  - atmBoundaryLayerInletK
  - atmBoundaryLayerInletEpsilon
  - atmBoundaryLayerInletOmega

## Benchmark (Physical phenomenon)

- The benchmark is an empty fetch computational
domain, steady-state RANS simulation.
- Flow characteristics:
  - External flow
  - The surface layer portion of the neutral-stratified
  equilibrium atmospheric boundary layer (no Ekman layer)
  - Dry air
  - Homogeneous, smooth terrain
  - Spatiotemporal-invariant aerodynamic roughness length
  - No displacement height
  - Newtonian, single-phase, incompressible, non-reacting
- Benchmark scenario:
  - Computational domain: (HW:Fig. 1)
  - Benchmark dataset: (HW:Fig. 6)
  (Obtained by the WebPlotDigitizer-4.2 (Rohatgi, 2019))

## Resources

    Computational study (tag:HW):
        Hargreaves, D. M., & Wright, N. G. (2007).
        On the use of the k–ε model in commercial CFD software
        to model the neutral atmospheric boundary layer.
        Journal of wind engineering and
        industrial aerodynamics, 95(5), 355-369.
        DOI:10.1016/j.jweia.2006.08.002

    Wind profile (tag:RQP):
        Richards, P. J., Quinn, A. D., & Parker, S. (2002).
        A 6 m cube in an atmospheric boundary layer flow-Part 2.
        Computational solutions. Wind and structures, 5(2_3_4), 177-192.
        DOI:10.12989/was.2002.5.2_3_4.177


# Numerics
## Physical modelling:

- The governing equations for:
  - Steady-state, Newtonian, single-phase, incompressible fluid flows,
    excluding any thermal chemical, electromagnetic and scalar interactions
- Mathematical approach for the turbulence modelling:
  - Reynolds-averaged Navier-Stokes simulation (RANS)
- Turbulence closure model:
  - kEpsilon and kOmegaSST linear eddy viscosity closure models
- The sets of input (HW:Table 1):
  - Reference height, Zref = 6 [m]
  - Aerodynamic roughness height, z0 = 0.01 [m]
  - Displacement height, d = 0 [m]
  - Reference mean wind speed, Uref = 10 [m/s]

## Computational domain modelling:

- Rectangular prism
- (x1, x2, x3) = (5000, 100, 500) [m] = (streamwise, spanwise, ground-normal) directions

## Computational domain discretisation:

- Spatial resolution:
  - (x1, x2, x3) = (500, 5, 50) [cells]
  - Refer to the `system/blockMeshDict` for the grading details
- Temporal resolution: Steady state

## Equation discretisation:

- Spatial derivatives and variables:
  - Convection: Second order
  - Others: Second order with various limiters
- Temporal derivatives and variables: First order

## Numerical boundary/initial conditions:

- Refer to `0.orig`

## Pressure-velocity coupling algorithm:

- SIMPLEC

## Linear solvers:

- Refer to `system/fvSolution`

## Initialisation and sampling:

- No initialisation/averaging
- Sampling at the end of the simulation via `system/sampleDict`
- Refer to `system/controlDict` for further details


<!------------------------------------------------------------------------- -->
