/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "PDRsetFields.H"
#include "PDRutilsInternal.H"
#include "mathematicalConstants.H"

#ifndef FULLDEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // A sign-corrected multiply
    // This is used for porosity of obstacle intersections
    inline static scalar COMBLK(const scalar a, const scalar b)
    {
        if (a < 0)
        {
            return -a * b;
        }

        return a * b;
    }


    // Obstacle satisfies some minimum size checks.
    // A volume check misses thin plates, so use area.
    // Thin sheet overlaps can be produced by touching objects
    // if the obs_extend parameter is > 0.
    inline static bool obsHasMinSize(const vector& span, const PDRparams& tol)
    {
        return
        (
            (cmptProduct(span) > tol.min_overlap_vol)
         &&
            (
                (span.x() * span.y() > tol.min_overlap_area)
             || (span.y() * span.z() > tol.min_overlap_area)
             || (span.z() * span.x() > tol.min_overlap_area)
            )
        );
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::PDRutils::one_d_overlap
(
    scalar xmin,
    scalar xmax,
    const PDRblock::location& grid,
    List<scalar>& olap,
    int *cmin, int *cmax,
    int *cfmin, int *cfmax
)
{
    // Looking at one coordinate direction, called x here, for something
    // that extends from xmin to xmax, calculate how much it overlaps
    // each cell in this direction. Result returned in 'olap' List is
    // the proportion of the grid step overlapped, i.e dimensionless.
    // First and last steps overlapped given by *cmin, *cmax
    // Ditto for shifted grid given by *cfmin, *cfmax.

    // Initially zero everywhere
    olap = Zero;

    if (olap.size() < grid.nPoints())
    {
        FatalErrorInFunction
            << "The overlap scratch array is too small, has "
            << olap.size() << " but needs " << grid.nPoints() << nl
            << exit(FatalError);
    }


    // No intersection with the box
    if (xmax <= grid.first() || grid.last() <= xmin)
    {
        // Mark as bad range, cannot iterate
        *cmin = 0;
        *cmax = -1;

        // Another bad range (cannot iterate) but for extra safety ensure
        // that (cfmin -> cmin) and (cmax -> cfmax) cannot iterate either
        *cfmin = 1;
        *cfmax = -2;
        return;
    }

    // Ensure search is within the (point) bounds
    xmin = grid.clip(xmin);
    xmax = grid.clip(xmax);

    // The begin/end of the obstacle
    *cmin = grid.findCell(xmin);
    *cmax = grid.findCell(xmax);

    for (label ix = *cmin; ix <= *cmax; ++ix)
    {
        olap[ix] = 1.0;
    }

    // Fixup ends
    if (*cmax == *cmin)
    {
        olap[*cmax] = (xmax - xmin) / grid.width(*cmax);
    }
    else
    {
        if (grid[*cmin] < xmin)
        {
            olap[*cmin] = (grid[*cmin+1] - xmin) / grid.width(*cmin);
        }

        if (xmax < grid[*cmax+1])
        {
            olap[*cmax] = (xmax - grid[*cmax]) / grid.width(*cmax);
        }
    }
    assert(olap[*cmax] >= 0.0);


    // Is xmin below/above the cell-centre (for virtual staggered-grid) ?
    *cfmin =
    (
        xmin < grid.C(*cmin)
      ? *cmin
      : Foam::min(*cmin+1, grid.nCells()-1)
    );

    // Is xmax below/above the cell-centre (for virtual staggered-grid) ?
    *cfmax =
    (
        xmax < grid.C(*cmax)
      ? *cmax
      : Foam::min(*cmax+1, grid.nCells()-1)
    );
}


/**************************************************************************************************/

void Foam::PDRutils::two_d_overlap
(
    const UList<scalar>& a_olap, label amin, label amax,
    const UList<scalar>& b_olap, label bmin, label bmax,
    SquareMatrix<scalar>& ab_olap
)
{
    // We go one over the relevant min/max limits since these values might be
    // used. If not, they would have been zeroed in one_d_overlap

    amin = Foam::max(0, amin-1);
    bmin = Foam::max(0, bmin-1);
    amax = Foam::min(a_olap.size()-1, amax+1);
    bmax = Foam::min(b_olap.size()-1, bmax+1);

    for (label ia = amin; ia <= amax; ++ia)
    {
        for (label ib = bmin; ib <= bmax; ++ib)
        {
            ab_olap(ia,ib) = a_olap[ia] * b_olap[ib];
        }
    }
}


/**************************************************************************************************/

void Foam::PDRutils::circle_overlap
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
)
{
    /* This routine calculates the proportion of each (two-dimensional) grid cell
       overlapped by the circle or angled rectangle. Coordinates are labelled a and b.
       On entry:
              ac, bc                coordinates of centre of circle or rectangle
              dia                   diameter of circle (zero for rectangle)
              theta, wa, wb parameters for rectangle
              agrid[]               locations of grid lines of a-grid
              amin, amax    first and last cells in a-grid overlapped by object
              (similarly for b)
            On exit:
              abolap                2-D array of (proportionate) area blockage by grid cell
              a_lblock              2-D array of (proportionate) blockage to a-direction flow
                                            (This will be area blockage when extruded in the third coordinate).
              a_count               (2-D array)The contribution of this object to the count of obstacles blocking
                                            a-direction flow. This is only non-zero if the object is inside the
                                            lateral boundaries of the cell. It is large negative if the cell is
                                            totally blocked in this direction.
              (similarly for b)
              c_drag                2-D array of tensor that will give tensor drag in each cell (when multiplied
                                            Cd, cylinder length, and 0.5 rho*U^2) Dimension: L.

       Note that this routine does not zero array elements outside the amin to amax, bmin to bmax area.
    */

    scalar count, a_lblk, b_lblk, perim, dummy;

    symmTensor2D vdrag(Zero);

    // Prevent stepping outside of the array when the obstacle is on the
    // upper boundary

    // Upper limit of inclusive range is nCells-1
    amin = Foam::max(0, amin);
    bmin = Foam::max(0, bmin);
    amax = Foam::min(amax, agrid.nCells()-1);
    bmax = Foam::min(bmax, bgrid.nCells()-1);

    for (label ia = amin; ia <= amax; ++ia)
    {
        // Cell-centred grid
        const scalar a1 = agrid[ia];
        const scalar a2 = agrid[ia+1];

        // Left-shifted staggered face grid (-1 addressing is OK)
        const scalar af1 = agrid.C(ia-1);
        const scalar af2 = agrid.C(ia);

        for (label ib = bmin; ib <= bmax; ++ib)
        {
            // Cell-centred grid
            const scalar b1 = bgrid[ib];
            const scalar b2 = bgrid[ib+1];

            // Left-shifted staggered face grid (-1 addressing is OK)
            const scalar bf1 = bgrid.C(ib-1);
            const scalar bf2 = bgrid.C(ib);

            // Do the centred cell
            if ( dia > 0.0 )
            {
                ab_olap(ia,ib) = inters_cy
                (
                    ac, bc, 0.5*dia, a1, a2, b1, b2, &perim,
                    &dummy, &dummy, &b_lblk, &a_lblk
                );
/* The last two arguments of the above call appear to be reversed, but the inters_cy routine returns
   the amount of overlap in the a and b direcvtions, which are the blockage to the b and a directions. */

/* abolap * cell area is area of cylinder in this cell. Divide by PI%D^2/4 to get proportion of cylinder in cell
   For whole cylinger c_drag should be = D, so multiply by D.                                                   */

                c_drag(ia,ib).xx() = c_drag(ia,ib).yy() = 4.0 * ab_olap(ia,ib) * (a2 - a1) * (b2 - b1) / dia / mathematical::pi;
                c_drag(ia,ib).xy() = Zero;
                c_count(ia,ib) = perim / (mathematical::pi * dia);

//******?????
                scalar area = (a2 - a1) * (b2 - b1);
                scalar rat = dia * dia / area - 1.5;
                if (rat > 0.0)
                {
                    scalar da = ac - 0.5 * (a1 + a2);
                    scalar db = bc - 0.5 * (b1 + b2);
                    scalar dc = std::hypot(da, db);
                    scalar rat1 = min(max((dc / sqrt(area) - 0.3) * 1.4, 0), 1);
                    scalar drg0 = c_drag(ia,ib).xx();
                    scalar drg1 = c_drag(ia,ib).yy();
                    scalar drg = std::hypot(drg0, drg1);
                    c_drag(ia,ib).xx() = drg * ( 1.0 - rat1 ) + drg * da*da/dc/dc * rat1;
                    c_drag(ia,ib).yy() = drg * ( 1.0 - rat1 ) + drg * db*db/dc/dc * rat1;
                    c_drag(ia,ib).xy() = drg * da*db/dc/dc *rat1;
                }
            }
            else
            {
                ab_olap(ia,ib) = inters_db( ac, bc, theta, wa, wb, a1, a2, b1, b2, &count, c_drag(ia,ib), &perim, &a_lblk, &b_lblk, &dummy, &dummy );
                c_count(ia,ib) = perim / ( wa + wb ) * 0.5;
            }
            ac_lblock(ia,ib) = a_lblk;
            bc_lblock(ia,ib) = b_lblk;
            ab_perim(ia,ib) = perim;

            // Do the a-shifted cell
            if ( dia > 0.0 ) // I.e. a cylinder, not a d.b.
            {
                if (ac >= af1 && ac < af2)
                {
                    // Only want to block one layer of faces
                    a_lblock(ia,ib) = l_blockage
                    (
                        ac, bc, 0.5*dia,
                        af1, af2, b1, b2, &count, &dummy, &dummy
                    );
                }
                inters_cy
                (
                    ac, bc, 0.5*dia,
                    af1, af2, b1, b2,
                    &perim, &count, &dummy, &dummy, &dummy
                );
            }
            else
            {
                inters_db
                (
                    ac, bc, theta, wa, wb, af1, af2, b1, b2,
                    &count, vdrag, &dummy, &a_lblk, &b_lblk, &dummy, &dummy
                );
                a_lblock(ia,ib) = a_lblk;
            }

            // Do the b-shifted cell
            if ( dia > 0.0 )
            {
                if (bc >= bf1 && bc < bf2)
                {
                    // Only want to block one layer of faces
                    b_lblock(ia,ib) = l_blockage
                    (
                        bc, ac, 0.5*dia, bf1, bf2, a1, a2,
                        &count, &(vdrag.yy()), &dummy
                    );
                }

                inters_cy
                (
                    ac, bc, 0.5*dia,
                    a1, a2, bf1, bf2,
                    &perim, &dummy, &count, &dummy, &dummy
                );
            }
            else
            {
                inters_db
                (
                    ac, bc, theta, wa, wb,
                    a1, a2, bf1, bf2,
                    &count, vdrag, &dummy, &a_lblk, &b_lblk, &dummy, &dummy
                );
                b_lblock(ia,ib) =  b_lblk;
            }
        }
    }

}  // End circle_overlap


/**************************************************************************************************/

scalar block_overlap
(
    DynamicList<PDRobstacle>& blocks,
    const labelRange& range,
    const scalar multiplier
)
{
    // Size information
    const label nBlock = range.size();

    // The return value
    scalar totVolume = 0;

    if (nBlock < 2) return 0;


    // Sort blocks by their x-position (with sortBias)
    labelList blkOrder;
    sortedOrder(blocks.slice(range), blkOrder);

    DynamicList<PDRobstacle> newBlocks;

    // Work through the sorted blocks
    for (label i1 = 0; i1 < nBlock-1; ++i1)
    {
        const PDRobstacle& blk1 = blocks[range[blkOrder[i1]]];

        // Upper coordinates
        const vector max1 = blk1.pt + blk1.span;

        // For second block start with the next one on the list, and
        // stop when we find the first one whose biased x-position
        // is beyond the end of the block1

        for (label i2 = i1 + 1; i2 < nBlock; ++i2)
        {
            const PDRobstacle& blk2 = blocks[range[blkOrder[i2]]];

            // Upper coordinates
            const vector max2 = blk2.pt + blk2.span;

            if (max1.x() <= blk2.x())
            {
                break;
            }

            if
            (
                max1.y() <= blk2.y()
             || max1.z() <= blk2.z()
             || max2.y() <= blk1.y()
             || max2.z() <= blk1.z()
             || (blk1.vbkge * blk2.vbkge <= 0)
            )
            {
                continue;
            }


            {
                PDRobstacle over;

                over.pt = max(blk1.pt, blk2.pt);
                over.span = min(max1, max2) - over.pt;

                assert(cmptProduct(over.span) > 0.0);

                // This routine should only have been called for all +ve o r all -ve obstacles
                assert(blk1.vbkge * blk2.vbkge > 0);
                /* At the first level of intersection, we create an obstacle of blockage -1 (if both objects solid)
                 to cancel out the double counting. (multiplier is 1).
                 ?? COMBLK does a (sign corrected) multiply; is this corrrect for porous obstacles?
                 Depends on how blockages were summed in the first place. In fact this -ve obstacle
                 concept only works if the blockages are summed??*/
                over.vbkge = - COMBLK( blk1.vbkge, blk2.vbkge ) * multiplier;
                over.xbkge = - COMBLK( blk1.xbkge, blk2.xbkge ) * multiplier;
                over.ybkge = - COMBLK( blk1.ybkge, blk2.ybkge ) * multiplier;
                over.zbkge = - COMBLK( blk1.zbkge, blk2.zbkge ) * multiplier;
                over.typeId = 81 + int(15 * multiplier);  // Not subsequently used

                if (obsHasMinSize(over.span, pars))
                {
                    // Obstacle satisfies some minimum size checks
                    totVolume -= over.volume();

                    newBlocks.append(over);
                }
            }
        }
    }

    blocks.append(std::move(newBlocks));

    return totVolume;
}


/**************************************************************************************************/

using namespace Foam::PDRutils;

scalar block_cylinder_overlap
(
    DynamicList<PDRobstacle>& blocks,
    const labelRange& range,
    const UList<PDRobstacle>& cylinders
)
{
    // Size information
    const label nBlock = range.size();
    const label nCyl = cylinders.size();

    // The return value
    scalar totVolume = 0;

    if (!nBlock || !nCyl) return 0;

    scalar area, a_lblk, b_lblk, dummy, a_centre, b_centre;
    symmTensor2D dum2;


    // Sort blocks and cylinders by their x-position (with sortBias)
    labelList blkOrder;
    sortedOrder(blocks.slice(range), blkOrder);

    labelList cylOrder;
    sortedOrder(cylinders, cylOrder);

    DynamicList<PDRobstacle> newBlocks;

    // Work through the sorted blocks
    for (label i1 = 0; i1 < nBlock; i1++)
    {
        const PDRobstacle& blk1 = blocks[range[blkOrder[i1]]];

        // Upper coordinates
        const vector max1 = blk1.pt + blk1.span;

        // Cyls whose end is before start of this block no longer
        // need to be considered

        label i2 = 0;
        while (i2 < nCyl-1 && cylinders[cylOrder[i2]] < blk1)
        {
            ++i2;
        }

        for (/*nil*/; i2 < nCyl; ++i2)
        {
            const PDRobstacle& cyl2 = cylinders[cylOrder[i2]];

            // Calculate overlap in axis direction; if zero continue.
            // Calculate 2-d overlap and c 0f g; if area zero continue.

            PDRobstacle over;


            switch (cyl2.orient)
            {
                case vector::Z:
                {
                    const scalar zm2 = cyl2.z() + cyl2.len();
                    if (blk1.z() > zm2 || cyl2.z() > max1.z()) continue;

                    if ( cyl2.dia() == 0.0 )
                    {
                        area = inters_db
                        (
                            cyl2.x(), cyl2.y(), cyl2.theta(), cyl2.wa, cyl2.wb,
                            blk1.x(), max1.x(),
                            blk1.y(), max1.y(),
                            &dummy, dum2, &dummy, &a_lblk, &b_lblk,
                            &a_centre, &b_centre
                        );
                    }
                    else
                    {
                        area = inters_cy
                        (
                            cyl2.x(), cyl2.y(), 0.5*cyl2.dia(),
                            blk1.x(), max1.x(),
                            blk1.y(), max1.y(),
                            &dummy, &dummy, &dummy, &dummy, &dummy
                        );
                        b_lblk = l_blockage
                        (
                            cyl2.x(), cyl2.y(), 0.5*cyl2.dia(),
                            blk1.x(), max1.x(),
                            blk1.y(), max1.y(),
                            &dummy, &dummy, &b_centre
                        );
                        a_lblk = l_blockage
                        (
                            cyl2.y(), cyl2.x(), 0.5*cyl2.dia(),
                            blk1.y(), max1.y(),
                            blk1.x(), max1.x(),
                            &dummy, &dummy, &a_centre
                        );
                    }
                    if (equal(area, 0)) continue;
                    assert(a_lblk >0.0);
                    assert(b_lblk >0.0);

                    // The intersection between a circle and a rectangle  can be an odd shape.
                    // We have its area. a_lblk and b_lblk are dimensions of enclosing rectangle
                    // and a_centre and b_centre its centre. We scale this rectangle down to
                    // the correct areacorrect area, as a rectangular approximation to the intersection.
                    const scalar ratio = std::sqrt( area / a_lblk / b_lblk );

                    a_lblk *= blk1.span.x() * ratio;
                    b_lblk *= blk1.span.y() * ratio;
                    assert(b_lblk >0.0);
                    assert(a_lblk >0.0);

                    over.x() = a_centre - 0.5 * a_lblk;
                    over.y() = b_centre - 0.5 * b_lblk;
                    over.z() = max(blk1.z(), cyl2.z());

                    over.span.x() = a_lblk;
                    over.span.y() = b_lblk;
                    over.span.z() = min(max1.z(), cyl2.z() + cyl2.len()) - over.z();
                    assert(over.x() > -200.0);
                    assert(over.x() < 2000.0);
                }
                break;

                case vector::Y:
                {
                    const scalar ym2 = cyl2.y() + cyl2.len();
                    if (blk1.y() > ym2 || cyl2.y() > max1.y()) continue;

                    if ( cyl2.dia() == 0.0 )
                    {
                        area = inters_db
                        (
                            cyl2.z(), cyl2.x(), cyl2.theta(), cyl2.wa, cyl2.wb,
                            blk1.z(), max1.z(),
                            blk1.x(), max1.x(),
                            &dummy, dum2, &dummy, &a_lblk, &b_lblk,
                            &a_centre, &b_centre
                        );
                    }
                    else
                    {
                        area = inters_cy
                        (
                            cyl2.z(), cyl2.x(), 0.5*cyl2.dia(),
                            blk1.z(), max1.z(),
                            blk1.x(), max1.x(),
                            &dummy, &dummy, &dummy, &dummy, &dummy
                        );

                        b_lblk = l_blockage
                        (
                            cyl2.z(), cyl2.x(), 0.5*cyl2.dia(),
                            blk1.z(), max1.z(),
                            blk1.x(), max1.x(),
                            &dummy, &dummy, &b_centre
                        );

                        a_lblk = l_blockage
                        (
                            cyl2.x(), cyl2.z(), 0.5*cyl2.dia(),
                            blk1.x(), max1.x(),
                            blk1.z(), max1.z(),
                            &dummy, &dummy, &a_centre
                        );
                    }

                    if (equal(area, 0)) continue;
                    assert(a_lblk >0.0);
                    assert(b_lblk >0.0);

                    // a_lblk and b_lblk are dimensions of enclosing rectangle.
                    // Need to scale to correct area
                    const scalar ratio = std::sqrt( area / a_lblk / b_lblk );
                    a_lblk *= blk1.span.z() * ratio;
                    b_lblk *= blk1.span.x() * ratio;

                    over.z() = a_centre - a_lblk * 0.5;
                    over.x() = b_centre - b_lblk * 0.5;
                    over.y() = max(blk1.y(), cyl2.y());

                    over.span.z() = a_lblk;
                    over.span.x() = b_lblk;
                    over.span.y() = min(max1.y(), cyl2.y() + cyl2.len()) - over.y();
                }
                break;

                case vector::X:
                {
                    const scalar xm2 = cyl2.x() + cyl2.len();
                    if (blk1.x() > xm2 || cyl2.x() > max1.x()) continue;

                    if ( cyl2.dia() == 0.0 )
                    {
                        area = inters_db
                        (
                            cyl2.y(), cyl2.z(), cyl2.theta(), cyl2.wa, cyl2.wb,
                            blk1.y(), max1.y(),
                            blk1.z(), max1.z(),
                            &dummy, dum2, &dummy, &a_lblk, &b_lblk,
                            &a_centre, &b_centre
                        );
                    }
                    else
                    {
                        area = inters_cy
                        (
                            cyl2.y(), cyl2.z(), 0.5*cyl2.dia(),
                            blk1.y(), max1.y(),
                            blk1.z(), max1.z(),
                            &dummy, &dummy, &dummy, &dummy, &dummy
                        );

                        b_lblk = l_blockage
                        (
                            cyl2.y(), cyl2.z(), 0.5*cyl2.dia(),
                            blk1.y(), max1.y(),
                            blk1.z(), max1.z(),
                            &dummy, &dummy, &b_centre
                        );

                        a_lblk = l_blockage
                        (
                            cyl2.z(), cyl2.y(), 0.5*cyl2.dia(),
                            blk1.z(), max1.z(),
                            blk1.y(), max1.y(),
                            &dummy, &dummy, &a_centre
                        );

                    }

                    if (equal(area, 0)) continue;
                    assert(a_lblk >0.0);
                    assert(b_lblk >0.0);

                    // a_lblk and b_lblk are dimensions of enclosing rectangle.
                    // Need to scale to correct area
                    const scalar ratio = std::sqrt( area / a_lblk / b_lblk );
                    assert(ratio >-10000.0);
                    assert(ratio <10000.0);
                    a_lblk *= blk1.span.y() * ratio;
                    b_lblk *= blk1.span.z() * ratio;

                    over.y() = a_centre - a_lblk * 0.5;
                    over.z() = b_centre - b_lblk * 0.5;
                    over.x() = max(blk1.x(), cyl2.x());

                    over.span.y() = a_lblk;
                    over.span.z() = b_lblk;
                    over.span.x() = min(max1.x(), cyl2.x() + cyl2.len()) - over.x();
                }
                break;
            }
            over.vbkge = over.xbkge = over.ybkge = over.zbkge = -1.0;
            over.typeId = PDRobstacle::IGNORE;

            assert(cmptProduct(over.span) > 0.0);
            assert(b_lblk >0.0);
            assert(a_lblk >0.0);
            assert(over.x() > -10000.0);

            if (obsHasMinSize(over.span, pars))
            {
                // Obstacle satisfies some minimum size checks
                totVolume -= over.volume();

                newBlocks.append(over);
            }
        }
    }

    blocks.append(std::move(newBlocks));

    return totVolume;
}


// ************************************************************************* //
