/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Preparation of fields for PDRFoam

\*---------------------------------------------------------------------------*/

#ifndef PDRsetFields_H
#define PDRsetFields_H

#include "PDRarrays.H"
#include "PDRblock.H"
#include "PDRmeshArrays.H"
#include "PDRobstacle.H"
#include "PDRpatchDef.H"
#include "PDRparams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


//YCYCLIC is set to 1 for running the test cases with a cyclic boundry condition in the y direction
//TWO_D is set to 1 for running the 2D test cases (no z direction) - usually same case as YCYCLIC
//Now specified in CAD_PDRDict and read in as globals.

// The program also labels intersection obstacles as types 96 and 86, but not then use these values
// Obstacles in group definitions have a multiple of 100 added to the type number

#define floatSMALL  1.e-10

// Default initial values of field variables, used outside congested area,
//and everywhere for uniform fields. They atre strings because same routines
//are used to create b.c.s for scalars and tensors.
#define DEFAULT_K       0.00015
#define DEFAULT_EPS     1e-5
#define DEFAULT_T       300
#define DEFAULT_P       100000
#define DEFAULT_SU      0.5
#define DEFAULT_LOBS    0.1   // Does not matter what it is outside congestion
                              // but zero would cause problems with Ep
#define DEFAULT_EP      0.01  // Gives length scale 0.1, calc. as (Xp-0.999)/Ep with Xp=1

// Boundary conditions on walls for all variables where it is not "zero_gradient"
#define K_WALL_FN       "kqRWallFunction"
#define EPS_WALL_FN     "epsilonWallFunction"
#define ALPHAT_WALL     "nutkWallFunction"
#define MUT_WALL_FN     "mutkWallFunction"
#define NUT_WALL_FN     "nutkWallFunction"

#define K_WALL_FN_LEGACY        "compressible::kqRWallFunction"
#define EPS_WALL_FN_LEGACY      "compressible::epsilonWallFunction"
#define ALPHAT_WALL_FN_LEGACY   "alphatWallFunction;\n\t\tPrt\t0.85"


// The following parameters are used to decide when there arMAX_Ne sufficient (parts of)
// obstacles ina cell for them to define the length scale of the generated turbulence.
#define MIN_AB_FOR_SIZE 0.002
#define MAX_VB_FOR_SIZE 0.9
#define COUNT_FOR_SIZE  0.1
#define MIN_COUNT_FOR_SIZE 0.05

// These define how blocked a face or cell has to be for removal from the mesh
//#define BLOCKED_CELL_PAR  0.05  //<-  Now pars.blockedCellPoros
//#define BLOCKED_FACE_PAR  0.95  //<-  Now pars.blockedFacePar


//- Calculate block/block overlaps
//
// Binary self-intersections are to be checked for blocks.
// Resulting negative blocks are appended to blocks.
// These new blocks have the opposite sign from input blocks, and
// blockage multiplied by multiplier.
//
// If the number of newly generated blocks is required, check the size
// of blocks on output vs input to see how many have been added.
//
// \param[in,out] blocks
// \param[in] range - the range within blocks to be examined
//
// \return overlap volume
scalar block_overlap
(
    DynamicList<PDRobstacle>& blocks,
    const labelRange& range,
    const scalar multiplier = 1.0
);


//- Calculate block/cylinder overlaps
//
// Binary intersections are to be checked for blocks and cylinders.
// Resulting negative blocks are appended to blocks.
// These new blocks have the opposite sign from input blocks, and
// blockage multiplied by multiplier.
//
// If the number of newly generated blocks is required, check the size
// of blocks on output vs input to see how many have been added.
//
// \param[in,out] arrp
// \param[in,out] blocks
// \param[in] range - the range within blocks to be examined
// \param[in] cylinders - the cylinders to be examined
//
//  \return overlap volume
scalar block_cylinder_overlap
(
    DynamicList<PDRobstacle>& blocks,
    const labelRange& range,
    const UList<PDRobstacle>& cylinders
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
