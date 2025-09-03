/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Namespace
    Foam::PDRutils

Description
    Utilities for PDR (eg, for setFields). Internal usage only.

    The C lineage of the original code is still evident in the use of
    pointers instead of references.
    This will be addressed in later versions of the code (2019-12).

SourceFiles
    PDRUtils.C

\*---------------------------------------------------------------------------*/

#ifndef PDRutilsInternal_H
#define PDRutilsInternal_H

#include "PDRutils.H"
#include "PDRarrays.H"
#include "PDRblock.H"
#include "symmTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRutils
{

//- Determine 1-D overlap locations for a geometric entity
//
// \param[in] xmin - min position of the geometric entity
// \param[in] xmax - max position of the geometric entity
// \param[in] grid - grid point information
// \param[out] olap - Fraction of cell-width with overlap
//     0 for no overlap, 1 for full overlap.
// \param[out] cmin - first cell index (inclusive) with overlap,
//     values in the range \c [0,nCells]
// \param[out] cmax - last cell index (inclusive) with overlap,
//     values in the range \c [0,nCells]
// \param[out] cfmin - first cell index (inclusive) with staggered face,
//     values in the range \c [0,nCells]
// \param[out] cfmax - last cell index (inclusive) with staggered face,
//     values in the range \c [0,nCells]
void one_d_overlap
(
    scalar xmin,
    scalar xmax,
    const PDRblock::location& grid,
    List<scalar>& olap,
    int *cmin, int *cmax,
    int *cfmin, int *cfmax
);


//- Combine two 1D overlaps.
//  Multiplying the two 1-d overlaps yields the proportion of each (2D) cell
//  that is covered.
//
//  \note We go one over the relevant min/max limits since these
//  values might be used.
//  The 1D arrays will have bee initially zeroed throughout.
void two_d_overlap
(
    const UList<scalar>& a_olap, label amin, label amax,
    const UList<scalar>& b_olap, label bmin, label bmax,
    SquareMatrix<scalar>& ab_olap
);


//- Calculate the proportion of each (two-dimensional) grid cell
//- overlapped by the circle or angled rectangle.
//
//  Coordinates are labelled a and b.
//
// \param[in] ac, bc  coordinates of centre of circle or rectangle
// \param[in] dia      diameter of circle (zero for rectangle)
// \param[in] theta, wa, wb     parameters for rectangle
// \param[in] amin, amax first and last cells in a-grid overlapped by object
// \param[in] agrid               locations of grid lines of a-grid
// \param[in] amin, amax first and last cells in b-grid overlapped by object
// \param[in] bgrid               locations of grid lines of b-grid
//
// \param[out] abolap
//     2-D array of (proportionate) area blockage by grid cell
// \param[out] a_lblock
//     2-D array of (proportionate) blockage to a-direction flow
//     (This will be area blockage when extruded in the third
//     coordinate).
//
// \param[out] a_count
//   2-D array The contribution of this object to the count of
// obstacles blocking a-direction flow. This is only non-zero if the
// object is inside the lateral boundaries of the cell. It is large
// negative if the cell is totally blocked in this direction.
//
//
// \param[out] c_drag
//
// 2-D array of tensor that will give tensor drag in each cell (when
// multiplied Cd, cylinder length, and 0.5 rho*U^2) Dimension: L.
//
// \note this routine does not zero array elements outside the amin
//  to amax, bmin to bmax area.
void circle_overlap
(
    scalar ac, scalar bc, scalar dia,
    scalar theta, scalar wa, scalar wb,
    const PDRblock::location& agrid, label amin, label amax,
    const PDRblock::location& bgrid, label bmin, label bmax,
    SquareMatrix<scalar>& ab_olap,
    SquareMatrix<scalar>& ab_perim,
    SquareMatrix<scalar>& a_lblock,
    SquareMatrix<scalar>& ac_lblock,
    SquareMatrix<scalar>& c_count,
    SquareMatrix<symmTensor2D>& c_drag,
    SquareMatrix<scalar>& b_lblock,
    SquareMatrix<scalar>& bc_lblock
);


//- Area of intersection between circle and rectangle.
//
// Calculates the area of intersection between the circle, centre (xc, yc), radius rad,
// and the rectangle with sides at x = x1 & x2, and y = y1 and y2.
//
// The return value is the fraction of the rectangle's area covered by the circle.
double inters_cy
(
    double xc,   //!< circle centre (x)
    double yc,   //!< circle centre (y)
    double rad,  //!< circle radius
    double x1, double x2,
    double y1, double y2,
    scalar* perim_p,
    scalar* x_proj_edge_p,
    scalar* y_proj_edge_p,
    scalar* x_overlap_p,
    scalar* y_overlap_p
);


//- The area overlap in the plane of a diagonal block and a cell.
//
// Calculates the overlap, in the plane of a diagonal block and a cell,
// plus blockage and drag parameters.
// Note that x and y herein may be any two of the three coordinates - would have been better not to label them x and y.
//
// On entry:
//   xc, yc                       Coordinates of axis of d.b.
//   theta, wa, wb        Angle and widths
//
//   The returned parameters will be multiplied by the length of the obstacle's intersection with
//   the third dimension of the 3-D cell to give this obstacle's contribution to the count, drag
//   and area blockages.
//   The return value is the area of intersection, which will multiply to volume blockage.
//
double inters_db
(
    double xc, double yc, double theta,
    double wa, double wb,
    double x1, double x2,
    double y1, double y2,
    scalar* count_p,
    symmTensor2D& vdrag, scalar* perim_p,
    scalar* x_lblk, scalar* y_lblk,
    scalar* x_centre_p, scalar* y_centre_p
);


/* Calculates the blockage to x-direction flow presented by the specified circle on
 the specified rectangle.
 Becomes the area blockage when extruded to in the third dimension.
 In other words, it is the projection on the y axis of the intersection between the
 circle and the rectangle.
 Returns fraction blocked
 Note that x and y in this routine may in fact be any two of the three dimensions.
 */
double l_blockage
(
    double xc, double yc, double rad,
    double x1, double x2,
    double y1, double y2,
    scalar* count_p, scalar* drag_p, scalar* centre_p
);


} // End namespace PDRutils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
