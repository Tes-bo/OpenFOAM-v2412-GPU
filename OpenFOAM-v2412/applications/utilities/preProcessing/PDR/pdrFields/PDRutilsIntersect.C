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

// Calculate the area of the sector of a circle whose ends are at
// (dxa, dya) and (dxb, dyb) relative to the centre. radsqu is radius
// squared.
//
// (We trust that this is consistent with the other parameters..)
inline static scalar sector
(
    scalar dxa, scalar dya,
    scalar dxb, scalar dyb
)
{
    scalar angle = (::atan2(dyb, dxb) - ::atan2(dya, dxa));

    if (angle < -1.0E-10)
    {
        angle += mathematical::twoPi;
    }

    return angle;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::PDRutils::inters_cy
(
    double xc, double yc, double rad,
    double x1, double x2,
    double y1, double y2,
    scalar* perim_p,
    scalar* x_proj_edge_p, scalar* y_proj_edge_p,
    scalar* x_overlap_p, scalar* y_overlap_p
)
{
    double angle, area, del;
    double x_int[6][2], y_int[6][2]; // Coordinates of intersections between the circle and sides of the rectangle.
    double x_arc[6][2], y_arc[6][2]; // Coordinates of end orc (within cell) for each quadrant
    double dx[6], dy[6];

    double x_olap_min = GREAT;
    double x_olap_max = -GREAT;
    double y_olap_min = GREAT;
    double y_olap_max = -GREAT;
    int n_vert, n_oppv;
    int no_intersection;

    const double dx1 = (x1 - xc);
    const double dx2 = (x2 - xc);
    const double dy1 = (y1 - yc);
    const double dy2 = (y2 - yc);

    const double dx1squ = dx1 * dx1;
    const double dx2squ = dx2 * dx2;
    const double dy1squ = dy1 * dy1;
    const double dy2squ = dy2 * dy2;
    const double radsqu = rad * rad;


    /* Going clockwise from (x1, y1), vertices are labelled 1,2,3,4, with 0 the same as 4
     and 5 the same as 1 (so that we can find the vertices on either side of any of 1 to 4).*/

    dx[1] = dx1; dy[1] = dy1;
    dx[2] = dx1; dy[2] = dy2;
    dx[3] = dx2; dy[3] = dy2;
    dx[4] = dx2; dy[4] = dy1;
    dx[0] = dx2; dy[0] = dy1;
    dx[5] = dx1; dy[5] = dy1;

    // The positions of the ends of the arcs, if these points are
    // inside the cell, they will be changed, if necessary, below.

    x_arc[2][0] = x_arc[1][1] = -rad; y_arc[2][0] = y_arc[1][1] =  0.0;
    x_arc[3][0] = x_arc[2][1] =  0.0; y_arc[3][0] = y_arc[2][1] =  rad;
    x_arc[4][0] = x_arc[3][1] =  rad; y_arc[4][0] = y_arc[3][1] =  0.0;
    x_arc[1][0] = x_arc[4][1] =  0.0; y_arc[1][0] = y_arc[4][1] = -rad;

    // We catch arcs that are entirely inside the rectangle
    // Note: this is wrong for a circle completely outside, but that
    // will be dealt with separately

    int arc_in[6] = { /* zero-initialied */ };
    arc_in[1] = (dx1 < -rad && dy1 < -rad) ? 1 : 0;
    arc_in[2] = (dx1 < -rad && dy2 >  rad) ? 1 : 0;
    arc_in[3] = (dx2 >  rad && dy2 >  rad) ? 1 : 0;
    arc_in[4] = (dx2 >  rad && dy1 < -rad) ? 1 : 0;

    // Work out which vertices are in the circle

    int vert_in[6];
    vert_in[1] = (dx1squ + dy1squ <= radsqu);
    vert_in[2] = (dx1squ + dy2squ <= radsqu);
    vert_in[3] = (dx2squ + dy2squ <= radsqu);
    vert_in[4] = (dx2squ + dy1squ <= radsqu);
    vert_in[0] = vert_in[4];
    vert_in[5] = vert_in[1];

    int n_in = 0;
    if (vert_in[1]) ++n_in;
    if (vert_in[2]) ++n_in;
    if (vert_in[3]) ++n_in;
    if (vert_in[4]) ++n_in;


    /* We now calculate the points of intersection of the circle with, successively,
     x=x1, y=y2, x=x2. y=y1.

     Where there are two intersections with one side, need to be careful about
     the order of the two points (i.e. clockwise round the rectangle) so that
     later on we get the right sector (short or long way round the circumference) */

    int n_int[6] = { /* zero-initialied */ };
    n_int[1] = 0;
    if ( dx1squ <= radsqu)
    {
        del = std::sqrt( radsqu - dx1squ);
        if ( ( ( -del ) <= dy2 ) && ( del >= dy1 ) )
        {
            x_int[1][0] = x_int[1][1] = dx1;
            if ( (-del ) > dy1 )
            {
                y_int[1][0] = -del;
                n_int[1]++;
                // This intersection will be an end of the 3rd- or 4th-quadrant arc
                if ( dx1 > 0.0 ) { x_arc[4][1] = dx1; y_arc[4][1] = -del; arc_in[4] = 1; }
                else             { x_arc[1][1] = dx1; y_arc[1][1] = -del; arc_in[1] = 1; }
            }
            if ( ( del ) < dy2 )
            {
                y_int[1][n_int[1]] =  del;
                n_int[1]++;
                if ( dx1 > 0.0 ) { x_arc[3][0] = dx1; y_arc[3][0] = del; arc_in[3] = 1; }
                else             { x_arc[2][0] = dx1; y_arc[2][0] = del; arc_in[2] = 1; }
            }
        }
    }

    n_int[2] = 0;
    if ( dy2squ <= radsqu)
    {
        del = std::sqrt( radsqu - dy2squ);
        if ( ( ( -del ) <= dx2 ) && ( del >= dx1 ) )
        {
            y_int[2][0] = y_int[2][1] = dy2;
            if ( (-del ) > dx1 )
            {
                x_int[2][0] = -del;
                n_int[2]++;
                if ( dy2 > 0.0 ) { x_arc[2][1] = -del; y_arc[2][1] = dy2; arc_in[2] = 1; }
                else             { x_arc[1][1] = -del; y_arc[1][1] = dy2; arc_in[1] = 1; }
            }
            if ( ( del ) < dx2 )
            {
                x_int[2][n_int[2]] =  del;
                n_int[2]++;
                if ( dy2 > 0.0 ) { x_arc[3][0] = del; y_arc[3][0] = dy2; arc_in[3] = 1; }
                else             { x_arc[4][0] = del; y_arc[4][0] = dy2; arc_in[4] = 1; }
            }
        }
    }

    n_int[3] = 0;
    if ( dx2squ <= radsqu)
    {
        del = std::sqrt( radsqu - dx2squ);
        if ( ( ( -del ) <= dy2 ) && ( del >= dy1 ) )
        {
            x_int[3][0] = x_int[3][1] = dx2;
            if ( ( del ) < dy2 )
            {
                y_int[3][0] = del;
                n_int[3]++;
                if ( dx2 > 0.0 ) { x_arc[3][1] = dx2; y_arc[3][1] = del; arc_in[3] = 1; }
                else             { x_arc[2][1] = dx2; y_arc[2][1] = del; arc_in[2] = 1; }
            }
            if ( (-del ) > dy1 )
            {
                y_int[3][n_int[3]] = -del;
                n_int[3]++;
                if ( dx2 > 0.0 ) { x_arc[4][0] = dx2; y_arc[4][0] = -del; arc_in[4] = 1; }
                else             { x_arc[1][0] = dx2; y_arc[1][0] = -del; arc_in[1] = 1; }
            }
        }
    }

    n_int[4] = 0;
    if ( dy1squ <= radsqu)
    {
        del = std::sqrt( radsqu - dy1squ);
        if ( ( ( -del ) <= dx2 ) && ( del >= dx1 ) )
        {
            y_int[4][0] = y_int[4][1] = dy1;
            if ( ( del ) < dx2 )
            {
                x_int[4][0] = del;
                n_int[4]++;
                if ( dy1 > 0.0 ) { x_arc[3][1] = del; y_arc[3][1] = dy1; arc_in[3] = 1; }
                else             { x_arc[4][1] = del; y_arc[4][1] = dy1; arc_in[4] = 1; }
            }
            if ( (-del ) > dx1 )
            {
                x_int[4][n_int[4]] = -del;
                n_int[4]++;
                if ( dy1 > 0.0 ) { x_arc[2][0] = -del; y_arc[2][0] = dy1; arc_in[2] = 1; }
                else             { x_arc[1][0] = -del; y_arc[1][0] = dy1; arc_in[1] = 1; }
            }
        }
    }

    n_int[0] = n_int[4];
    n_int[5] = n_int[1];

    y_int[0][0] = y_int[0][1] = dy1;
    x_int[0][0] = x_int[4][0];
    x_int[0][1] = x_int[4][1];
    x_int[5][0] = x_int[5][1] = dx1;
    y_int[5][0] = y_int[1][0];
    y_int[5][1] = y_int[1][1];

    /* There are five separate cases, depending of the number of vertices inside the circle */
    switch ( n_in )
    {
        case 0:
        {
            /* We start with the whole area of the circle, and then subtract any bits that stick out. */
            area =  radsqu * mathematical::pi;
            *perim_p = mathematical::twoPi * rad;
            no_intersection = true;
            for (n_vert = 1; n_vert < 5; n_vert++)
            {
                assert(n_int[n_vert] != 1);
                if (n_int[n_vert] == 2)
                {
                    /* The area of the bit to be subtracted is a sector minus a triangle. */
                    no_intersection = false;
                    angle = sector( x_int[n_vert][1], y_int[n_vert][1], x_int[n_vert][0], y_int[n_vert][0]);
                    area -= angle * radsqu * 0.5;
                    *perim_p -= angle * rad;
                    /* Two trinagles specified here, but one has zero area. */
                    area += ( - ( y_int[n_vert][1] - y_int[n_vert][0] ) * x_int[n_vert][0]
                              + ( x_int[n_vert][1] - x_int[n_vert][0] ) * y_int[n_vert][0] ) / 2.0;
                }
            }
            /* Need to allow for when the circle is completely out side the rectanglle
             by checking if the centre is outside the rectangle                                           */
            if ( no_intersection )
            {
                if ( (dx1>0) ||(dx2<0) || (dy1>0) || (dy2<0) )
                {
                    *perim_p = *x_proj_edge_p = *y_proj_edge_p = 0.0;
                    area = *x_overlap_p = *y_overlap_p = 0.0;
                    return area;
                }
            }

            break;
        }

        case 1:
        {
            /* Find which vertex is inside */
            n_vert = 1;
            while ( !vert_in[n_vert] ) { n_vert++; assert( n_vert < 5 ); }
            assert( n_int[n_vert-1] == 1 );
            if ( n_int[n_vert] != 1 )
            {
                assert( n_int[n_vert] == 1 );
            }
            angle = sector( x_int[n_vert-1][0], y_int[n_vert-1][0], x_int[n_vert][0], y_int[n_vert][0]);
            area = angle * radsqu * 0.5;
            *perim_p = angle * rad;
            /* We subtract (or add) two triangles; the other two evaluate to zero */
            area -= ( - ( x_int[n_vert][0]   - dx[n_vert] ) * dy[n_vert]
                      + ( x_int[n_vert-1][0] - dx[n_vert] ) * dy[n_vert]
                      + ( y_int[n_vert][0]   - dy[n_vert] ) * dx[n_vert]
                      - ( y_int[n_vert-1][0] - dy[n_vert] ) * dx[n_vert] ) / 2.0;

            break;
        }

        case 2:
        {
            /* This time n_vert is the number of the side which is completely inside the circle */
            n_vert = 1;
            while ( !(vert_in[n_vert] && vert_in[n_vert+1]) ) { n_vert++; assert( n_vert < 5 ); }
            assert( n_int[n_vert-1] == 1 );
            assert( n_int[n_vert+1] == 1 );
            angle = sector( x_int[n_vert-1][0], y_int[n_vert-1][0], x_int[n_vert+1][0], y_int[n_vert+1][0]);
            area = angle * radsqu * 0.5;
            *perim_p = angle * rad;
            /* We subtract (or add) three triangles; the other three evaluate to zero */
            area += ( ( x_int[n_vert+1][0] - dx[n_vert+1] ) * dy[n_vert+1]
                      - ( x_int[n_vert-1][0] - dx[n_vert] ) * dy[n_vert]
                      - ( y_int[n_vert+1][0] - dy[n_vert+1] ) * dx[n_vert+1]
                      + ( y_int[n_vert-1][0] - dy[n_vert] ) * dx[n_vert]
                      + ( dx[n_vert+1] -dx[n_vert] ) * dy[n_vert]
                      - ( dy[n_vert+1] -dy[n_vert] ) * dx[n_vert] ) / 2.0;

            switch ( n_vert )
            {
                case 1: x_olap_min = dx1; break;
                case 2: y_olap_max = dy2; break;
                case 3: x_olap_max = dx2; break;
                case 4: y_olap_min = dy1; break;
            }

            break;
        }

        case 3:
        {
            /* Find which vertex is NOT inside */
            n_vert = 1;
            while ( vert_in[n_vert] ) { n_vert++; assert( n_vert < 5 ); }
            assert( n_int[n_vert-1] == 1 );
            assert( n_int[n_vert] == 1 );
            n_oppv = (n_vert + 2) % 4;
            angle = sector( x_int[n_vert][0], y_int[n_vert][0], x_int[n_vert-1][0], y_int[n_vert-1][0]);
            area = angle * radsqu * 0.5;
            *perim_p = angle * rad;
            /* We subtract (or add) four triangles; the other four evaluate to zero */
            area += ( - ( x_int[n_vert][0]   - dx[n_vert+1] ) * dy[n_vert+1]
                      + ( x_int[n_vert-1][0] - dx[n_vert-1] ) * dy[n_vert-1]
                      + ( y_int[n_vert][0]   - dy[n_vert+1] ) * dx[n_vert+1]
                      - ( y_int[n_vert-1][0] - dy[n_vert-1] ) * dx[n_vert-1]
                      + ( dx[n_oppv] -dx[n_vert+1] ) * dy[n_oppv]
                      - ( dx[n_oppv] -dx[n_vert-1] ) * dy[n_oppv]
                      - ( dy[n_oppv] -dy[n_vert+1] ) * dx[n_oppv]
                      + ( dy[n_oppv] -dy[n_vert-1] ) * dx[n_oppv] ) / 2.0;

            x_olap_min = dx1;
            y_olap_max = dy2;
            x_olap_max = dx2;
            y_olap_min = dy1;

            break;
        }

        case 4:
        {
            /* Easy! We have the whole rectangle.  */
            area = *x_overlap_p = *y_overlap_p = 1.0; // Normalised
            *perim_p = *x_proj_edge_p = *y_proj_edge_p = 0.0;
            return area;

            break;
        }
    }

    // The area may be very small negative by rounding errors
    assert(area >=-1.0E-4);
    if (area < 0.0) area = 0.0;
    /* Return the overlap as a fraction of the rectangle's area. */
    area /= ( (x2 - x1 ) * ( y2 - y1 ) );

    // Sum the parts of the circumference that are inside the circle, projected onto axes
    *x_proj_edge_p =
    (
        (y_arc[1][1] - y_arc[1][0]) * arc_in[1]
      + (y_arc[2][1] - y_arc[2][0]) * arc_in[2]
      + (y_arc[3][0] - y_arc[3][1]) * arc_in[3]
      + (y_arc[4][0] - y_arc[4][1]) * arc_in[4]
    );

    *y_proj_edge_p =
    (
        (x_arc[1][0] - x_arc[1][1]) * arc_in[1]
      + (x_arc[2][1] - x_arc[2][0]) * arc_in[2]
      + (x_arc[3][1] - x_arc[3][0]) * arc_in[3]
      + (x_arc[4][0] - x_arc[4][1]) * arc_in[4]
    );

    if (arc_in[1])
    {
        x_olap_min = min(x_olap_min, x_arc[1][1]);
        x_olap_max = max(x_olap_max, x_arc[1][0]);
        y_olap_min = min(y_olap_min, y_arc[1][0]);
        y_olap_max = max(y_olap_max, y_arc[1][1]);
    }
    if (arc_in[2])
    {
        x_olap_min = min(x_olap_min, x_arc[2][0]);
        x_olap_max = max(x_olap_max, x_arc[2][1]);
        y_olap_min = min(y_olap_min, y_arc[2][0]);
        y_olap_max = max(y_olap_max, y_arc[2][1]);
    }
    if (arc_in[3])
    {
        x_olap_min = min(x_olap_min, x_arc[3][0]);
        x_olap_max = max(x_olap_max, x_arc[3][1]);
        y_olap_min = min(y_olap_min, y_arc[3][1]);
        y_olap_max = max(y_olap_max, y_arc[3][0]);
    }
    if (arc_in[4])
    {
        x_olap_min = min(x_olap_min, x_arc[4][1]);
        x_olap_max = max(x_olap_max, x_arc[4][0]);
        y_olap_min = min(y_olap_min, y_arc[4][1]);
        y_olap_max = max(y_olap_max, y_arc[4][0]);
    }

    *x_overlap_p = ( x_olap_max - x_olap_min ) / ( x2 - x1 );
    *y_overlap_p = ( y_olap_max - y_olap_min ) / ( y2 - y1 );
    assert ( *x_overlap_p >= -floatSMALL );
    assert ( *y_overlap_p >= -floatSMALL );

    return area;
} // End intersect


// ************************************************************************* //

double Foam::PDRutils::l_blockage
(
    double xc, double yc, double rad,
    double x1, double x2,
    double y1, double y2,
    scalar* count_p, scalar* drag_p, scalar* centre_p
)
{
    double xi = 0.0, lb, lb1, lb2, del;
    bool within = true; // Indicates that the the intersection does not overlap the ends of the line

    /* xi is the side we need to calc. intersections with */
    if      ( xc < x1 ) { xi = x1; }
    else if ( xc > x2 ) { xi = x2; }

    if ( xi == 0.0 )
    {
        del = rad;     // The relevant lowest ( or highest) point is at end of vertical radius
    }
    else                // The relevant lowest ( or highest) point at intersection with x = xi
    {
        del = rad*rad - ( xi - xc ) * ( xi - xc );
        if ( del < 0.0 ) { del = 0.0; }                  // No intersection
        else             { del = std::sqrt(del); }
    }

    if ( ( yc + del ) > y2 ) { lb2  = y2; within = false; } else { lb2  = yc + del; }
    if ( ( yc - del ) < y1 ) { lb1 = y1; within = false; } else { lb1 = yc - del; }

    lb = (lb2 - lb1) / (y2 - y1);
    *centre_p = (lb2 + lb1) * 0.5;

    if ( lb < 0.0 )  lb = 0.0;

    /* *count_p is 0 if the circle overlaps either y-side of the rectangle,
     1 if the circle is entirely inside the rectangle
     reduced if it overlaps x-sides.
     A negative value indicates total blockage*/
    if ( within && (lb > 0.0) )
    {
        *count_p = 1.0;
        if ( ( xc - rad ) < x1 ) *count_p -= 0.5;
        if ( ( xc + rad ) > x2 ) *count_p -= 0.5;
    }
    else
    {
        *count_p = 0.0;
    }
    *drag_p = lb * 1.2; //*drag_p = lb * CD_ROUND;
    if ( lb > 0.99 ) { *count_p = -1000.0; *drag_p = 1000.0; }
    assert(lb >-100.0);
    return lb;
}// End l_blockage


// ************************************************************************* //

double Foam::PDRutils::inters_db
(
    double xc, double yc, double theta,
    double wa, double wb,
    double x1, double x2,
    double y1, double y2,
    scalar* count_p,
    symmTensor2D& vdrag, scalar* perim_p,
    scalar* x_lblk_p, scalar* y_lblk_p,
    scalar* x_centre_p, scalar* y_centre_p
)
{
    double  x_int[6][2], y_int[6][2]; // Coordinates of intersections between the circle and sides of the rectangle.
    double  area, lpa, lpb, len;

    double  m = ::tan( theta );
    double  cth = ::cos( theta );
    double  sth = ::sin( theta );

    double  was  = wa * sth * 0.5;
    double  wac  = wa * cth * 0.5;
    double  wbs  = wb * sth * 0.5;
    double  wbc  = wb * cth * 0.5;
    double  waos = wa / sth * 0.5;
    double  waoc = wa / cth * 0.5;
    double  wbos = wb / sth * 0.5;
    double  wboc = wb / cth * 0.5;

    double xb[6], yb[6], xp1, xp2, yp1, yp2;

    double  dx1 = (x1 - xc);
    double  dx2 = (x2 - xc);
    double  dy1 = (y1 - yc);
    double  dy2 = (y2 - yc);

    *count_p = 0;

// The vertices of the rectangle (all coordinates relative to centre of rectangle)
    xb[1] = -wac - wbs;
    yb[1] = -was + wbc;
    xb[3] =  wac + wbs;
    yb[3] =  was - wbc;
    xb[2] =  wac - wbs;
    yb[2] =  was + wbc;
    xb[4] = -wac + wbs;
    yb[4] = -was - wbc;

    // First parameter of x_int or y_int determines which side of the cell we intersecting with
    //  Second parameter 0 is first intersection, 1 is second, going clockwise

    if ( xb[1] < dx1 )  // i.e. if left corner of block is to the left of x1
    {
        // Where one of lower sides of block intersects with x=x1
        // Innermost max determines which intersection is the genuine one
        // (not if whole block is to left of x1)
        y_int[1][0] = min(max(max(dx1 * m - wboc, -dx1 / m - waos), dy1), dy2);
        // Upper intersection
        y_int[1][1] = min(max(min(dx1 * m + wboc, -dx1 / m + waos), dy1), dy2);
    }
    else
    {
        y_int[1][1] = dy1;
        y_int[1][0] = dy2;
        // We add a quarter to count for each vertex inside the cell
        if ( (yb[1] > dy1) && (yb[1] < dy2) ) // ?? Seems inefficient ??
        { *count_p += 0.25; }
    }
    if ( xb[3] > dx2 )
    {
        y_int[3][1] = min(max(max(dx2 * m - wboc, -dx2 / m - waos), dy1), dy2);
        y_int[3][0] = min(max(min(dx2 * m + wboc, -dx2 / m + waos), dy1), dy2);
    }
    else
    {
        y_int[3][0] = dy1;
        y_int[3][1] = dy2;
        if (yb[3] > dy1 && yb[3] < dy2)
        {
            *count_p += 0.25;
        }
    }
    if (yb[2] > dy2)
    {
        x_int[2][0] = min(max(max(dy2 / m - wbos, -dy2 * m - waoc), dx1), dx2);
        x_int[2][1] = min(max(min(dy2 / m + wbos, -dy2 * m + waoc), dx1), dx2);
    }
    else
    {
        x_int[2][0] = dx2;
        x_int[2][1] = dx1;
        if ( (xb[2] > dx1) && (xb[2] < dx2) )
        { *count_p += 0.25; }
    }
    if ( yb[4] < dy1 )
    {
        x_int[4][1] = min(max(max(dy1 / m - wbos, -dy1 * m - waoc ), dx1), dx2);
        x_int[4][0] = min(max(min(dy1 / m + wbos, -dy1 * m + waoc ), dx1), dx2);
    }
    else
    {
        x_int[4][1] = dx2;
        x_int[4][0] = dx1;
        if ( (xb[4] > dx1) && (xb[4] < dx2) )
        { *count_p += 0.25; }
    }

    y_int[0][0] = y_int[0][1] = dy1;
    x_int[0][0] = x_int[4][0];
    x_int[0][1] = x_int[4][1];
    x_int[5][0] = x_int[5][1] = dx1;
    y_int[5][0] = y_int[1][0];
    y_int[5][1] = y_int[1][1];


// We can now define a smaller enclosing rectangle

    xp1 = min(x_int[2][0], x_int[4][1]);  // Leftmost of the intersections with top and bottom of cell
    if ( yb[1] > dy1 && yb[1] < dy2 )       xp1 = min(xp1, xb[1] ); // left corner of block
    xp1 = max(xp1, dx1);  // Make sure it is not to the left of x1

    yp2 = max(y_int[1][1], y_int[3][0] );
    if ( xb[2] > dx1 && xb[2] < dx2 )  yp2 = max(yp2, yb[2] );
    yp2 = min(yp2, dy2);

    xp2 = max(x_int[2][1], x_int[4][0] );
    if ( yb[3] > dy1 && yb[3] < dy2 )  xp2 = max(xp2, xb[3] );
    xp2 = min(xp2, dx2);

    yp1 = min(y_int[1][0], y_int[3][1]);
    if ( xb[4] > dx1 && xb[4] < dx2 )  yp1 = min(yp1, yb[4] );
    yp1 = max(yp1, dy1 );

    // Conveniently, the dimensions of the enclosing rectangle give us the line blockages
    *x_lblk_p = (xp2 - xp1 ) / (x2 - x1 );
    if ( *x_lblk_p < 0.0 ) { *x_lblk_p = 0.0; *count_p = 0.0; };  // ?? Better to trap no intersection earlier??
    *y_lblk_p = (yp2 - yp1 ) / (y2 - y1 );
    if ( *y_lblk_p < 0.0 ) { *y_lblk_p = 0.0; *count_p = 0.0; };

    *x_centre_p = xc + (xp2 + xp1 ) * 0.5;
    *y_centre_p = yc + (yp2 + yp1 ) * 0.5;

    *perim_p = lpa = lpb = 0.0;;
    area = (xp2 - xp1 ) * ( yp2 - yp1 );
    {
        double  dxx, dyy;
        // Lower left
        dyy = max(0.0, min(yb[1], y_int[1][0]) - yp1);
        dxx = min(xb[4], x_int[0][1] ) - xp1;

        if ( ( dxx * dyy) > 0.0 )
        {
            area -= dxx * dyy * 0.5;
            len = std::hypot(dxx, dyy);
            lpa += len * 0.5;
            *perim_p += len;
        }
        // Upper left
        dxx = max(0.0, min(xb[2], x_int[2][0]) - xp1);
        dyy = yp2 - max(yb[1], y_int[1][1] );
        if ( ( dxx * dyy) > 0.0 )
        {
            area -= dxx * dyy * 0.5;
            len = std::hypot(dxx, dyy);
            lpb += len * 0.5;
            *perim_p += len;
        }
        // Upper right
        dyy = max(0.0, yp2 - max(yb[3], y_int[3][0]));
        dxx = xp2 - max(xb[2], x_int[2][1] );
        if ( ( dxx * dyy) > 0.0 )
        {
            area -= dxx * dyy * 0.5;
            len = std::hypot(dxx, dyy);
            lpa += len * 0.5;
            *perim_p += len;
        }
        // Lower right
        dxx = max(0.0, xp2 - max(xb[4], x_int[4][0]));
        dyy = min(yb[3], y_int[3][1] ) - yp1;
        if ( ( dxx * dyy) > 0.0 )
        {
            area -= dxx * dyy * 0.5;
            len = std::hypot(dxx, dyy);
            lpb += len * 0.5;
            *perim_p += len;
        }

    }

    vdrag.xx() = lpa * cth * cth + lpb * sth * sth;
    vdrag.xy() = lpa * cth * sth - lpb * sth * cth;
    vdrag.yy() = lpa * sth * sth + lpb * cth * cth;

    return area / ( (x2 - x1 ) * ( y2 - y1 ) );
} // End inters_db


// ************************************************************************* //
