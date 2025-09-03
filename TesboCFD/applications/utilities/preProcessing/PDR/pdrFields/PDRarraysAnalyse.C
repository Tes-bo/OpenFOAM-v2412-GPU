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
#include "PDRobstacle.H"
#include "PDRpatchDef.H"
#include "PDRutils.H"
#include "PDRutilsInternal.H"
#include "ListOps.H"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#ifndef FULLDEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

using namespace Foam;
using namespace Foam::PDRutils;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Cell blockage.
//
// Use simple sum, because this will eventually be divided by a notional
// number of rows to give a per-row blockage ratio
//
//      b is the (fractional) area blockage. f is 1 or 0.5 to split between ends
//      Thus if b isclose to 1, the obstacle is totally blocking the cell in this direction,
//      and we could modify the behaviour if we wish.
inline static void add_blockage_c
(
    scalar& a,
    bool& blocked,
    const scalar b,
    const scalar f = 1.0
)
{
    a += b * f;
    if (b > pars.blockageNoCT)
    {
        blocked = true;
    }
}


// Face blockage
//
// Adds more area blockage to existing amount by assuming partial overlap,
// i.e. multiplying porosities.
//
// Simple addition if the existing amount is negative, because negative
// blocks (summed first) should just cancel out part of positive blocks.
inline static void add_blockage_f
(
    scalar& a,
    const scalar b,
    bool isHole
)
{
    if (a > 0.0)
    {
        // Both positive
        a = 1.0 - (1.0 - a) * (1.0 - b);
    }
    else if (b < pars.blockedFacePar || isHole)
    {
        // Add until it eventually becomes positive
        a += b;
    }
    else
    {
        // If one obstacle blocks face, face is blocked, regardless of
        // overlap calculations, unless an input negative obstacle makes a
        // hole in it
        a = b;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRarrays::addCylinder(const PDRobstacle& obs)
{
    if (equal(obs.vbkge, 0))
    {
        return;
    }

    if (isNull(block()))
    {
        FatalErrorInFunction
            << "No PDRblock set" << nl
            << exit(FatalError);
    }

    const PDRblock& pdrBlock = block();
    const PDRblock::location& xgrid = pdrBlock.grid().x();
    const PDRblock::location& ygrid = pdrBlock.grid().y();
    const PDRblock::location& zgrid = pdrBlock.grid().z();

    scalarList& xoverlap = overlap_1d.x();
    scalarList& yoverlap = overlap_1d.y();
    scalarList& zoverlap = overlap_1d.z();

    int     cxmin, cxmax, cymin, cymax, czmin, czmax;
    int     cfxmin, cfxmax, cfymin, cfymax, cfzmin, cfzmax;

    scalar rad_a, rad_b;
    vector area_block(Zero);

    if (obs.typeId == PDRobstacle::CYLINDER)
    {
        rad_a = rad_b = 0.5*obs.dia();
    }
    else
    {
        rad_a = 0.5*(obs.wa * ::cos(obs.theta()) + obs.wb * ::sin(obs.theta()));
        rad_b = 0.5*(obs.wa * ::sin(obs.theta()) + obs.wb * ::cos(obs.theta()));
    }

    switch (obs.orient)
    {
        case vector::Z:
        {
            // Determine the part of the grid (potentially) covered by this obstacle.
            one_d_overlap
            (
                obs.x() - rad_a, obs.x() + rad_a,
                pdrBlock.grid().x(),
                xoverlap, &cxmin, &cxmax, &cfxmin, &cfxmax
            ); assert(cxmax >=0);

            one_d_overlap
            (
                obs.y() - rad_b, obs.y() + rad_b,
                pdrBlock.grid().y(),
                yoverlap, &cymin, &cymax, &cfymin, &cfymax
            ); assert(cymax >=0);

            one_d_overlap
            (
                obs.z(), obs.z() + obs.len(),
                pdrBlock.grid().z(),
                zoverlap, &czmin, &czmax, &cfzmin, &cfzmax
            ); assert(czmax >=0);

            // The area of intersection with each 2D cell in an x-y plane.
            //  a corresponds to x, and b to y.
            circle_overlap
            (
                obs.x(), obs.y(), obs.dia(),
                obs.theta(), obs.wa, obs.wb,
                pdrBlock.grid().x(), cxmin, cfxmax,
                pdrBlock.grid().y(), cymin, cfymax,
                aboverlap, abperim, a_lblock, ac_lblock, c_count, c_drag,
                b_lblock, bc_lblock
            );


            for (label ix = cxmin; ix <= cfxmax; ix++)
            {
                for (label iy = cymin; iy <= cfymax; iy++)
                {
                    for (label iz = czmin; iz <= cfzmax; iz++)
                    {
                        const scalar vol_block = aboverlap(ix,iy) * zoverlap[iz];
                        v_block(ix,iy,iz) += vol_block;
                        surf(ix,iy,iz) += abperim(ix,iy) * zoverlap[iz] * zgrid.width(iz);

                        // In the 2D case, the ends of the cylinder appear to
                        // be in the cell even when not, making surf and
                        // obs_size wrong. So leave out ends.

                        if (!pars.two_d && (iz == czmin || iz == czmax))
                        {
                            // End cells
                            const scalar both_ends_fac = (czmin == czmax ? 2.0 : 1.0);

                            surf(ix,iy,iz) += aboverlap(ix,iy)
                                * xgrid.width(ix) * ygrid.width(iy) * both_ends_fac;
                        }

                        const scalar temp = c_count(ix,iy) * zoverlap[iz];

                        obs_count(ix,iy,iz) += temp;
                        sub_count(ix,iy,iz).z() += temp;

                        if (!pars.two_d && (iz == czmin || iz == czmax))
                        {
                            // End faces
                            const scalar both_ends_fac = (czmin == czmax ? 2.0 : 1.0);

                            sub_count(ix,iy,iz).z() -= temp * both_ends_fac / 2.0;
                        }

                        // Keep the blockage and drag of round obst separate
                        // from the sharp for the moment because different
                        // blockage ratio corrections will be needed later.
                        //
                        // Only the relevant diagonal elements of drag
                        // are stored; other are zero.

                        area_block.x() = ac_lblock(ix,iy) * zoverlap[iz];
                        area_block.y() = bc_lblock(ix,iy) * zoverlap[iz];

                        // Do not want blockage and drag across the end
                        // except at perimeter.
                        if (aboverlap(ix,iy) < pars.blockedFacePar)
                        {
                            if (obs.typeId == PDRobstacle::CYLINDER)
                            {
                                drag_r(ix,iy,iz).x() += c_drag(ix,iy).xx() * zoverlap[iz] / xgrid.width(ix) / ygrid.width(iy);
                                drag_r(ix,iy,iz).y() += c_drag(ix,iy).yy() * zoverlap[iz] / xgrid.width(ix) / ygrid.width(iy);

                                add_blockage_c(area_block_r(ix,iy,iz).x(), dirn_block(ix,iy,iz).x(), area_block.x(), 1.0);
                                add_blockage_c(area_block_r(ix,iy,iz).y(), dirn_block(ix,iy,iz).y(), area_block.y(), 1.0);
                            }
                            else
                            {
                                drag_s(ix,iy,iz).xx() += c_drag(ix,iy).xx() * zoverlap[iz] / xgrid.width(ix) / ygrid.width(iy);
                                drag_s(ix,iy,iz).yy() += c_drag(ix,iy).yy() * zoverlap[iz] / xgrid.width(ix) / ygrid.width(iy);

                                add_blockage_c(area_block_s(ix,iy,iz).x(), dirn_block(ix,iy,iz).x(), area_block.x(), 1.0);
                                add_blockage_c(area_block_s(ix,iy,iz).y(), dirn_block(ix,iy,iz).y(), area_block.y(), 1.0);
                            }
                        }
                        // Here we accumulate 1/betai - 1.
                        betai_inv1(ix,iy,iz).x() += vol_block / (1.0 - area_block.x() + floatSMALL);
                        betai_inv1(ix,iy,iz).y() += vol_block / (1.0 - area_block.y() + floatSMALL);
                        betai_inv1(ix,iy,iz).z() += vol_block / (1.0 - aboverlap(ix,iy) + floatSMALL);

                        // The off-diagonal elements of drag are stored in "drag" for BOTH round and sharp ostacles
                        drag_s(ix,iy,iz).xy() += c_drag(ix,iy).xy() * zoverlap[iz] / xgrid.width(ix) / ygrid.width(iy);

                        add_blockage_f(face_block(ix,iy,iz).x(), a_lblock(ix,iy) * zoverlap[iz], hole_in_face(ix,iy,iz).x());
                        add_blockage_f(face_block(ix,iy,iz).y(), b_lblock(ix,iy) * zoverlap[iz], hole_in_face(ix,iy,iz).y());
                    }

                    // Face blockage in the z direction
                    if (czmin == czmax)
                    {
                        // Does not span cell. Block nearest face.
                        add_blockage_f(face_block(ix,iy,cfzmax).z(), aboverlap(ix,iy), hole_in_face(ix,iy,cfzmax).z());
                    }
                    else
                    {
                        // In at least two cells.
                        // Block first and last faces overlapped
                        add_blockage_f(face_block(ix,iy,czmin+1).z(), aboverlap(ix,iy), hole_in_face(ix,iy,czmin+1).z());
                        if (czmax > czmin+1)
                        {
                            add_blockage_f(face_block(ix,iy,czmax).z(), aboverlap(ix,iy), hole_in_face(ix,iy,czmax).z());
                        }
                    }

                    // z_block is used to work out the blockage ratio for each
                    // "row" of sub-grid obstacles so this cylinder should
                    // not have any eeffct in cells that it completely spans.
                    // Hence statement commented out below and replaced by
                    // two lines after this loop. That longitudinal clockage
                    // goes into new array z_aling_block

                    for (label iz = czmin+1; iz < czmax; ++iz)
                    {
                        // Internal only
                        along_block(ix,iy,iz).z() += aboverlap(ix,iy);
                    }

                    // Longitudinal drag only applies at ends of cylinder.
                    // If cylinder spans more than 1 cell, apply half at each
                    // end.

                    drag_s(ix,iy,czmin).zz() += aboverlap(ix,iy) * xgrid.width(ix) * ygrid.width(iy) / 2.0;
                    drag_s(ix,iy,czmax).zz() += aboverlap(ix,iy) * xgrid.width(ix) * ygrid.width(iy) / 2.0;

                    add_blockage_c(area_block_s(ix,iy,czmin).z(), dirn_block(ix,iy,czmin).z(), aboverlap(ix,iy), 0.5);
                    add_blockage_c(area_block_s(ix,iy,czmax).z(), dirn_block(ix,iy,czmax).z(), aboverlap(ix,iy), 0.5);
                }
            }

            break;
        }

        case vector::X:  // orientation
        {
            // x-direction cylinder. a,b are y,z.
            one_d_overlap
            (
                obs.x(), obs.x() + obs.len(),
                pdrBlock.grid().x(),
                xoverlap, &cxmin, &cxmax, &cfxmin, &cfxmax
            ); assert(cxmax >=0);

            one_d_overlap
            (
                obs.y() - rad_a, obs.y() + rad_a,
                pdrBlock.grid().y(),
                yoverlap, &cymin, &cymax, &cfymin, &cfymax
            ); assert(cymax >=0);

            one_d_overlap
            (
                obs.z() - rad_b, obs.z() + rad_b,
                pdrBlock.grid().z(),
                zoverlap, &czmin, &czmax, &cfzmin, &cfzmax
            ); assert(czmax >=0);

            circle_overlap
            (
                obs.y(), obs.z(), obs.dia(),
                obs.theta(), obs.wa, obs.wb,
                pdrBlock.grid().y(), cymin, cfymax,
                pdrBlock.grid().z(), czmin, cfzmax,
                aboverlap, abperim, a_lblock, ac_lblock, c_count, c_drag,
                b_lblock, bc_lblock
            );


            for (label iy = cymin; iy <= cfymax; iy++)
            {
                for (label iz = czmin; iz <= cfzmax; iz++)
                {
                    for (label ix = cxmin; ix <= cxmax; ix++)
                    {
                        const scalar vol_block = aboverlap(iy,iz) * xoverlap[ix];
                        v_block(ix,iy,iz) += vol_block;
                        surf(ix,iy,iz) += abperim(iy,iz) * xoverlap[ix] * xgrid.width(ix);

                        if (ix == cxmin || ix == cxmax)
                        {
                            // End cells
                            const scalar both_ends_fac = (cxmin == cxmax ? 2.0 : 1.0);

                            surf(ix,iy,iz) += aboverlap(iy,iz)
                                * ygrid.width(iy) * zgrid.width(iz) * both_ends_fac;
                        }

                        const scalar temp = c_count(iy,iz) * xoverlap[ix];
                        obs_count(ix,iy,iz) += temp;
                        sub_count(ix,iy,iz).x() += temp;

                        if (ix == cfxmin || ix == cfxmax)
                        {
                            // End faces
                            const scalar both_ends_fac = (cfxmin == cfxmax ? 2.0 : 1.0);

                            sub_count(ix,iy,iz).x() -= temp * both_ends_fac / 2.0;
                        }

                        area_block.y() = ac_lblock(iy,iz) * xoverlap[ix];
                        area_block.z() = bc_lblock(iy,iz) * xoverlap[ix];

                        // Do not want blockage and drag across the end
                        // except at perimeter.
                        if (aboverlap(iy,iz) < pars.blockedFacePar)
                        {
                            if (obs.typeId == PDRobstacle::CYLINDER)
                            {
                                drag_r(ix,iy,iz).y() += c_drag(iy,iz).xx() * xoverlap[ix] / ygrid.width(iy) / zgrid.width(iz);
                                drag_r(ix,iy,iz).z() += c_drag(iy,iz).yy() * xoverlap[ix] / ygrid.width(iy) / zgrid.width(iz);

                                add_blockage_c(area_block_r(ix,iy,iz).y(), dirn_block(ix,iy,iz).y(), area_block.y(), 1.0);
                                add_blockage_c(area_block_r(ix,iy,iz).z(), dirn_block(ix,iy,iz).z(), area_block.z(), 1.0);
                            }
                            else
                            {
                                drag_s(ix,iy,iz).yy() += c_drag(iy,iz).xx() * xoverlap[ix] / ygrid.width(iy) / zgrid.width(iz);
                                drag_s(ix,iy,iz).zz() += c_drag(iy,iz).yy() * xoverlap[ix] / ygrid.width(iy) / zgrid.width(iz);

                                add_blockage_c(area_block_s(ix,iy,iz).y(), dirn_block(ix,iy,iz).y(), area_block.y(), 1.0);
                                add_blockage_c(area_block_s(ix,iy,iz).z(), dirn_block(ix,iy,iz).z(), area_block.z(), 1.0);
                            }
                        }
                        betai_inv1(ix,iy,iz).y() += vol_block / (1.0 - area_block.y() + floatSMALL);
                        betai_inv1(ix,iy,iz).z() += vol_block / (1.0 - area_block.z() + floatSMALL);
                        betai_inv1(ix,iy,iz).x() += vol_block / (1.0 - aboverlap(iy,iz) + floatSMALL);

                        // The off-diagonal elements of drag are stored in "drag" for BOTH round and sharp ostacles
                        drag_s(ix,iy,iz).yz() += c_drag(iy,iz).xy() * xoverlap[ix] / ygrid.width(iy) / zgrid.width(iz);

                        add_blockage_f(face_block(ix,iy,iz).y(), a_lblock(iy,iz) * xoverlap[ix], hole_in_face(ix,iy,iz).y());
                        add_blockage_f(face_block(ix,iy,iz).z(), b_lblock(iy,iz) * xoverlap[ix], hole_in_face(ix,iy,iz).z());
                    }
                    if (cxmin == cxmax)
                    {
                        // Does not span cell. Block nearest face.
                        add_blockage_f(face_block(cfxmax,iy,iz).x(), aboverlap(iy,iz), hole_in_face(cfxmax,iy,iz).x());
                    }
                    else
                    {
                        // In at least two cells.
                        // Block first and last faces overlapped
                        add_blockage_f(face_block(cxmin+1,iy,iz).x(), aboverlap(iy,iz), hole_in_face(cxmin+1,iy,iz).x());
                        if (cxmax > cxmin+1)
                        {
                            add_blockage_f(face_block(cxmax,iy,iz).x(), aboverlap(iy,iz), hole_in_face(cxmax,iy,iz).x());
                        }
                    }

                    for (label ix = cxmin+1; ix < cxmax; ++ix)
                    {
                        // Internal only
                        along_block(ix,iy,iz).x() += aboverlap(iy,iz);
                    }
                    drag_s(cxmin,iy,iz).xx() += aboverlap(iy,iz) * ygrid.width(iy) * zgrid.width(iz) / 2.0;
                    drag_s(cxmax,iy,iz).xx() += aboverlap(iy,iz) * ygrid.width(iy) * zgrid.width(iz) / 2.0;

                    add_blockage_c(area_block_s(cxmin,iy,iz).x(), dirn_block(cxmin,iy,iz).x(), aboverlap(iy,iz), 0.5);
                    add_blockage_c(area_block_s(cxmax,iy,iz).x(), dirn_block(cxmax,iy,iz).x(), aboverlap(iy,iz), 0.5);
                }
            }

            break;
        }

        case vector::Y:  // orientation
        {
            // y-direction cylinder. a,b are z,x.
            one_d_overlap
            (
                obs.x() - rad_b, obs.x() + rad_b,
                pdrBlock.grid().x(),
                xoverlap, &cxmin, &cxmax, &cfxmin, &cfxmax
            ); assert(cxmax >=0);

            one_d_overlap
            (
                obs.y(), obs.y() + obs.len(),
                pdrBlock.grid().y(),
                yoverlap, &cymin, &cymax, &cfymin, &cfymax
            ); assert(cymax >=0);

            one_d_overlap
            (
                obs.z() - rad_a, obs.z() + rad_a,
                pdrBlock.grid().z(),
                zoverlap, &czmin, &czmax, &cfzmin, &cfzmax
            ); assert(czmax >=0);

            circle_overlap
            (
                obs.z(), obs.x(), obs.dia(),
                obs.theta(), obs.wa, obs.wb,
                pdrBlock.grid().z(), czmin, cfzmax,
                pdrBlock.grid().x(), cxmin, cfxmax,
                aboverlap, abperim, a_lblock, ac_lblock, c_count, c_drag,
                b_lblock, bc_lblock
            );


            for (label iz = czmin; iz <= cfzmax; iz++)
            {
                for (label ix = cxmin; ix <= cfxmax; ix++)
                {
                    for (label iy = cymin; iy <= cymax; iy++)
                    {
                        const scalar vol_block = aboverlap(iz,ix) * yoverlap[iy];
                        v_block(ix,iy,iz) += vol_block;
                        surf(ix,iy,iz) += abperim(iz,ix) * yoverlap[iy] * ygrid.width(iy);

                        if (iy == cymin || iy == cymax)
                        {
                            // End cells
                            const scalar both_ends_fac = (cymin == cymax ? 2.0 : 1.0);

                            surf(ix,iy,iz) += aboverlap(iz,ix)
                                * zgrid.width(iz) * xgrid.width(ix) * both_ends_fac;
                        }

                        const scalar temp = c_count(iz,ix) * yoverlap[iy];

                        obs_count(ix,iy,iz) += temp;
                        sub_count(ix,iy,iz).y() += temp;

                        if (iy == cfymin || iy == cfymax)
                        {
                            // End faces
                            const scalar both_ends_fac = (cfymin == cfymax ? 2.0 : 1.0);

                            sub_count(ix,iy,iz).y() -= temp * both_ends_fac / 2.0;
                        }

                        area_block.z() = ac_lblock(iz,ix) * yoverlap[iy];
                        area_block.x() = bc_lblock(iz,ix) * yoverlap[iy];

                        // Do not want blockage and drag across the end
                        // except at perimeter.
                        if (aboverlap(iz,ix) < pars.blockedFacePar)
                        {
                            if (obs.typeId == PDRobstacle::CYLINDER)
                            {
                                drag_r(ix,iy,iz).z() += c_drag(iz,ix).xx() * yoverlap[iy] / zgrid.width(iz) / xgrid.width(ix);
                                drag_r(ix,iy,iz).x() += c_drag(iz,ix).yy() * yoverlap[iy] / zgrid.width(iz) / xgrid.width(ix);

                                add_blockage_c(area_block_r(ix,iy,iz).z(), dirn_block(ix,iy,iz).z(), area_block.z(), 1.0);
                                add_blockage_c(area_block_r(ix,iy,iz).x(), dirn_block(ix,iy,iz).x(), area_block.x(), 1.0);
                            }
                            else
                            {
                                drag_s(ix,iy,iz).zz() += c_drag(iz,ix).xx() * yoverlap[iy] / zgrid.width(iz) / xgrid.width(ix);
                                drag_s(ix,iy,iz).xx() += c_drag(iz,ix).yy() * yoverlap[iy] / zgrid.width(iz) / xgrid.width(ix);

                                add_blockage_c(area_block_s(ix,iy,iz).z(), dirn_block(ix,iy,iz).z(), area_block.z(), 1.0);
                                add_blockage_c(area_block_s(ix,iy,iz).x(), dirn_block(ix,iy,iz).x(), area_block.x(), 1.0);
                            }
                        }
                        betai_inv1(ix,iy,iz).z() += vol_block / (1.0 - area_block.z() + floatSMALL);
                        betai_inv1(ix,iy,iz).x() += vol_block / (1.0 - area_block.x() + floatSMALL);
                        betai_inv1(ix,iy,iz).y() += vol_block / (1.0 - aboverlap(iz,ix) + floatSMALL);

                        // The off-diagonal elements of drag are stored in "drag" for BOTH round and sharp obstacles
                        drag_s(ix,iy,iz).xz() += c_drag(iz,ix).xy() * yoverlap[iy] / zgrid.width(iz) / xgrid.width(ix);

                        add_blockage_f(face_block(ix,iy,iz).z(), a_lblock(iz,ix) * yoverlap[iy], hole_in_face(ix,iy,iz).z());
                        add_blockage_f(face_block(ix,iy,iz).x(), b_lblock(iz,ix) * yoverlap[iy], hole_in_face(ix,iy,iz).x());
                    }
                    if (cymin == cymax)
                    {
                        // Does not span cell. Block nearest face.
                        add_blockage_f(face_block(ix,cfymax,iz).y(), aboverlap(iz,ix), hole_in_face(ix,cfymax,iz).y());
                    }
                    else
                    {
                        // In at least two cells.
                        // Block first and last faces overlapped
                        add_blockage_f(face_block(ix,cymin+1,iz).y(), aboverlap(iz,ix), hole_in_face(ix,cymin+1,iz).y());
                        if (cymax > cymin+1)
                        {
                            add_blockage_f(face_block(ix,cymax,iz).y(), aboverlap(iz,ix), hole_in_face(ix,cymax,iz).y());
                        }
                    }

                    for (label iy = cymin+1; iy < cymax; ++iy)
                    {
                        // Internal only
                        along_block(ix,iy,iz).y() += aboverlap(iz,ix);
                    }

                    drag_s(ix,cymin,iz).yy() += aboverlap(iz,ix) * zgrid.width(iz) * xgrid.width(ix) / 2.0;
                    drag_s(ix,cymax,iz).yy() += aboverlap(iz,ix) * zgrid.width(iz) * xgrid.width(ix) / 2.0;

                    add_blockage_c(area_block_s(ix,cymin,iz).y(), dirn_block(ix,cymin,iz).y(), aboverlap(iz,ix), 0.5);
                    add_blockage_c(area_block_s(ix,cymax,iz).y(), dirn_block(ix,cymax,iz).y(), aboverlap(iz,ix), 0.5);
                }
            }

            break;
        }

        default:  // orientation
        {
            Info<< "Unexpected orientation " << int(obs.orient) << nl;
            break;
        }
    }
}


void Foam::PDRarrays::addBlockage
(
    const PDRobstacle& obs,
    DynamicList<PDRpatchDef>& patches,
    const int volumeSign
)
{
    // The volumeSign indicates whether this pass is for negative or positive
    // obstacles. Both if 0.
    if
    (
        equal(obs.vbkge, 0)
     || (volumeSign < 0 && obs.vbkge >= 0)
     || (volumeSign > 0 && obs.vbkge < 0)
    )
    {
        return;
    }

    if (isNull(block()))
    {
        FatalErrorInFunction
            << "No PDRblock set" << nl
            << exit(FatalError);
    }

    const PDRblock& pdrBlock = block();
    const PDRblock::location& xgrid = pdrBlock.grid().x();
    const PDRblock::location& ygrid = pdrBlock.grid().y();
    const PDRblock::location& zgrid = pdrBlock.grid().z();

    scalarList& xoverlap = overlap_1d.x();
    scalarList& yoverlap = overlap_1d.y();
    scalarList& zoverlap = overlap_1d.z();


    // 0 will be used later for faces found to be blocked.
    // 2 is used for wall-function faces.
    label patchNum = PDRpatchDef::LAST_PREDEFINED;

    // Only the part of the panel that covers full cell faces will be used
    // so later should keep the panel in the list plus a hole (-ve obstacle) for
    // part that will become blowoff b.c.

    int indir = 0;

    // Panel or patch
    const bool isPatch =
    (
        (obs.typeId == PDRobstacle::LOUVRE_BLOWOFF && obs.blowoff_type > 0)
     || (obs.typeId == PDRobstacle::RECT_PATCH)
    );

    if (isPatch)
    {
        const auto& identifier = obs.identifier;

        const auto spc = identifier.find_first_of(" \t\n\v\f\r");

        const word patchName = word::validate(identifier.substr(0, spc));

        patchNum = ListOps::find
        (
            patches,
            [=](const PDRpatchDef& p){ return patchName == p.patchName; },
            1  // skip 0 (blocked face)
        );

        if (patchNum < 1)
        {
            // The patch name was not already in the list
            patchNum = patches.size();

            patches.append(PDRpatchDef(patchName));
        }


        PDRpatchDef& p = patches[patchNum];

        if (obs.typeId == PDRobstacle::RECT_PATCH)
        {
            indir = obs.inlet_dirn;
            p.patchType = 0;
        }
        else
        {
            p.patchType = obs.blowoff_type;
            p.blowoffPress = obs.blowoff_press;
            p.blowoffTime = obs.blowoff_time;
            if (obs.span.x() < 0.01)
            {
                indir = 1;
            }
            else if (obs.span.y() < 0.01)
            {
                indir = 2;
            }
            else if (obs.span.z() < 0.01)
            {
                indir = 3;
            }
            else
            {
                FatalErrorInFunction
                    << "Blowoff panel should have zero thickness" << nl
                    << exit(FatalError);
            }
        }
    }

    int cxmin, cxmax, cfxmin, cfxmax;
    one_d_overlap
    (
        obs.x(), obs.x() + obs.span.x(),
        pdrBlock.grid().x(),
        xoverlap, &cxmin, &cxmax, &cfxmin, &cfxmax
    ); assert(cxmax >=0);

    int cymin, cymax, cfymin, cfymax;
    one_d_overlap
    (
        obs.y(), obs.y() + obs.span.y(),
        pdrBlock.grid().y(),
        yoverlap, &cymin, &cymax, &cfymin, &cfymax
    ); assert(cymax >=0);

    int czmin, czmax, cfzmin, cfzmax;
    one_d_overlap
    (
        obs.z(), obs.z() + obs.span.z(),
        pdrBlock.grid().z(),
        zoverlap, &czmin, &czmax, &cfzmin, &cfzmax
    ); assert(czmax >=0);

    two_d_overlap(xoverlap, cxmin, cxmax, yoverlap, cymin, cymax, aboverlap);


    const scalar vbkge = obs.vbkge;
    const scalar xbkge = obs.xbkge;
    const scalar ybkge = obs.ybkge;
    const scalar zbkge = obs.zbkge;

    // Compensate for double-counting of drag if two edges in same cell
    const vector double_f
    (
        ((cxmin == cxmax) ? 0.5 : 1.0),
        ((cymin == cymax) ? 0.5 : 1.0),
        ((czmin == czmax) ? 0.5 : 1.0)
    );

    for (label ix = cxmin; ix <= cfxmax; ix++)
    {
        const scalar xov = xoverlap[ix];

        scalar area, cell_area, temp;

        for (label iy = cymin; iy <= cfymax; iy++)
        {
            const scalar yov = yoverlap[iy];

            for (label iz = czmin; iz <= cfzmax; iz++)
            {
                const scalar zov = zoverlap[iz];

                if
                (
                    isPatch
                 &&
                    (
                        (indir == -1 && ix == cfxmin)
                     || (indir == 1  && ix == cfxmax)
                     || (indir == -2 && iy == cfymin)
                     || (indir == 2  && iy == cfymax)
                     || (indir == -3 && iz == cfzmin)
                     || (indir == 3  && iz == cfzmax)
                    )
                )
                {
                    /* Type RECT_PATCH (16) exists to set all faces it covers to be in a particular patch
                     usually an inlet or outlet.
                     ?? If the face not on a cell boundary, this will move it to the lower-cordinate
                     face of the relevant cell. It should be at the face of teh volume blocked by
                     the obstacle. But, if that is not at a cell boundary, the obstacle will be
                     putting blockage in front of teh vent, so we should be checking that it is
                     at a cell boundary. */

                    switch (indir) // Face orientation
                    {
                        // X
                        case -1:
                        case 1:
                            if (yov * zov > pars.blockedFacePar)
                            {
                                face_patch(ix,iy,iz).x() = patchNum;
                            }
                            break;

                        // Y
                        case -2:
                        case 2:
                            if (zov * xov > pars.blockedFacePar)
                            {
                                face_patch(ix,iy,iz).y() = patchNum;
                            }
                            break;

                        // Z
                        case -3:
                        case 3:
                            if (xov * yov > pars.blockedFacePar)
                            {
                                face_patch(ix,iy,iz).z() = patchNum;
                            }
                            break;
                    }
                }  // End of code for patch


                const scalar vol_block = aboverlap(ix,iy) * zov * vbkge;
                v_block(ix,iy,iz) += vol_block;

                // These are the face blockages
                if ((ix > cxmin && ix <= cxmax) || (cxmin == cxmax && ix == cfxmax))
                {
                    temp = yov * zov * xbkge;

                    // Has -ve volumeSign only when processing user-defined
                    // -ve blocks
                    if (volumeSign < 0)
                    {
                        hole_in_face(ix,iy,iz).x() = true;
                    }
                    add_blockage_f(face_block(ix,iy,iz).x(), temp, hole_in_face(ix,iy,iz).x());
                    if (temp > pars.blockedFacePar && ! hole_in_face(ix,iy,iz).x() && !isPatch)
                    {
                        // Put faces of block in wall patch
                        face_patch(ix,iy,iz).x() = PDRpatchDef::WALL_PATCH;
                    }
                }
                if ((iy > cymin && iy <= cymax) || (cymin == cymax && iy == cfymax))
                {
                    temp = zov * xov * ybkge;
                    if (volumeSign < 0)
                    {
                        hole_in_face(ix,iy,iz).y() = true;
                    }
                    add_blockage_f(face_block(ix,iy,iz).y(), temp, hole_in_face(ix,iy,iz).y());
                    if (temp > pars.blockedFacePar && ! hole_in_face(ix,iy,iz).y() && !isPatch)
                    {
                        face_patch(ix,iy,iz).y() = PDRpatchDef::WALL_PATCH;
                    }
                }
                if ((iz > czmin && iz <= czmax) || (czmin == czmax && iz == cfzmax))
                {
                    temp = xov * yov * zbkge;
                    if (volumeSign < 0)
                    {
                        hole_in_face(ix,iy,iz).z() = true;
                    }
                    add_blockage_f(face_block(ix,iy,iz).z(), temp, hole_in_face(ix,iy,iz).z());
                    if (temp > pars.blockedFacePar && ! hole_in_face(ix,iy,iz).z() && !isPatch)
                    {
                        face_patch(ix,iy,iz).z() = PDRpatchDef::WALL_PATCH;
                    }
                }

                // These are the interior blockages
                /* As for cylinders, blockage that extends longitudinally all the way through the cell
                 should not be in x_block etc., but it does go into new arrays along_block etc. */
                area = yov * zov * xbkge; // Note this is fraction of cell c-s area
                if (ix < cxmin || ix > cxmax)
                {}
                else if (ix > cxmin && ix < cxmax && xbkge >= 1.0)
                {
                    along_block(ix,iy,iz).x() += area;
                    betai_inv1(ix,iy,iz).x() += vol_block / (1.0 - area + floatSMALL);
                }
                else if (ix == cxmin || ix == cxmax)
                {
                    // If front and back of the obstacle are not in the
                    // same cell, put half in each
                    const scalar double_f = (cxmin == cxmax ? 1.0 : 0.5);

                    add_blockage_c(area_block_s(ix,iy,iz).x(), dirn_block(ix,iy,iz).x(), area, double_f);
                    cell_area = (ygrid.width(iy) * zgrid.width(iz));
                    surf(ix,iy,iz) += double_f * area * cell_area;
                    betai_inv1(ix,iy,iz).x() += vol_block / (1.0 - area + floatSMALL);

                    // To get Lobs right for a grating, allow for the fact that it is series of bars by having additional count
                    if (obs.typeId == PDRobstacle::GRATING && obs.orient == vector::X)
                    {
                        //  * cell_area to make dimensional, then / std::sqrt(cell_area) to get width
                        temp = area * std::sqrt(cell_area) / obs.slat_width - 1.0;
                        if (temp > 0.0) { grating_count(ix,iy,iz).x() += temp; }
                    }
                }

                /* As for cylinders, blockage that extends longitudinally all the way through the cell
                 should not be in x_block etc., but it does go into new arrays along_block etc. */
                area = zov * xov * ybkge;
                if (iy < cymin || iy > cymax)
                {}
                else if (iy > cymin && iy < cymax && ybkge >= 1.0)
                {
                    along_block(ix,iy,iz).y() += area;
                    betai_inv1(ix,iy,iz).y() += vol_block / (1.0 - area + floatSMALL);
                }
                else if (iy == cymin || iy == cymax)
                {
                    // If front and back of the obstacle are not in the
                    // same cell, put half in each
                    const scalar double_f = (cymin == cymax ? 1.0 : 0.5);

                    add_blockage_c(area_block_s(ix,iy,iz).y(), dirn_block(ix,iy,iz).y(), area, double_f);
                    cell_area = (zgrid.width(iz) * xgrid.width(ix));
                    surf(ix,iy,iz) += double_f * area * cell_area;
                    betai_inv1(ix,iy,iz).y() += vol_block / (1.0 - area + floatSMALL);

                    if (obs.typeId == PDRobstacle::GRATING && obs.orient == vector::Y)
                    {
                        //  * cell_area to make dimensional, then / std::sqrt(cell_area) to get width
                        temp = area * std::sqrt(cell_area) / obs.slat_width - 1.0;
                        if (temp > 0.0) { grating_count(ix,iy,iz).y() += temp; }
                    }
                }

                area = xov * yov * zbkge;
                if (iz < czmin || iz > czmax)
                {}
                else if (iz > czmin && iz < czmax && zbkge >= 1.0)
                {
                    along_block(ix,iy,iz).z() += area;
                    betai_inv1(ix,iy,iz).z() += vol_block / (1.0 - area + floatSMALL);
                }
                else if (iz == czmin || iz == czmax)
                {
                    // If front and back of the obstacle are not in the
                    // same cell, put half in each
                    const scalar double_f = (czmin == czmax ? 1.0 : 0.5);

                    add_blockage_c(area_block_s(ix,iy,iz).z(), dirn_block(ix,iy,iz).z(), area, double_f);
                    cell_area = (xgrid.width(ix) * ygrid.width(iy));
                    surf(ix,iy,iz) += double_f * area * cell_area;
                    betai_inv1(ix,iy,iz).z() += vol_block / (1.0 - area + floatSMALL);

                    if (obs.typeId == PDRobstacle::GRATING && obs.orient == vector::Z)
                    {
                        //  * cell_area to make dimensional, then / std::sqrt(cell_area) to get width
                        temp = area * std::sqrt(cell_area) / obs.slat_width - 1.0;
                        if (temp > 0.0) { grating_count(ix,iy,iz).z() += temp; }
                    }
                }
            }
        }
    }

    if (obs.typeId == PDRobstacle::RECT_PATCH)
    {
        // Was only needed to set face_patch values
        return;
    }


    /*  A narrow obstacle completely crossing the cell adds 1 to the count for the transverse directions
    If all four edges are not in the relevant cell, we take 1/4 for each edge that is in.
    If it does not totally cross the cell, then the count is reduced proportionately
    If it is porous, then the count is reduced proportionately.
    ?? Should it be? At least this is consistent with effect going smoothly to zero as porosity approaches 1 ??

    Note that more than 1 can be added for one obstacle, if sides and ends are in the cell.

    We do all the x-aligned edges first, both for y-blockage and z-blockage.                                */

    /* Intersection obstacles can have more edges than the intersecting obstacles that
     generated them, so not good to incorporate these in N and drag.  */

    for (label ix = cxmin; ix <= cxmax; ix++)
    {
        // Factor of 0.25 because 4 edges to be done
        scalar olap25 = 0.25 * xoverlap[ix];

        const scalar temp =
            ((pars.noIntersectN && vbkge < 0.334) ? 0 : (olap25 * vbkge));

        obs_count(ix,cymin,czmin) += temp;
        obs_count(ix,cymax,czmin) += temp;
        obs_count(ix,cymin,czmax) += temp;
        obs_count(ix,cymax,czmax) += temp;

        sub_count(ix,cymin,czmin).x() += temp;
        sub_count(ix,cymax,czmin).x() += temp;
        sub_count(ix,cymin,czmax).x() += temp;
        sub_count(ix,cymax,czmax).x() += temp;

        // The 0.25 becomes 0.5 to allow for front/back faces in drag direction
        olap25 *= 2.0;

        drag_s(ix,cymin,czmin).yy() += zoverlap[czmin] * double_f.z() * olap25 * ybkge / ygrid.width(cymin);
        drag_s(ix,cymax,czmin).yy() += zoverlap[czmin] * double_f.z() * olap25 * ybkge / ygrid.width(cymax);
        drag_s(ix,cymin,czmax).yy() += zoverlap[czmax] * double_f.z() * olap25 * ybkge / ygrid.width(cymin);
        drag_s(ix,cymax,czmax).yy() += zoverlap[czmax] * double_f.z() * olap25 * ybkge / ygrid.width(cymax);

        drag_s(ix,cymin,czmin).zz() += yoverlap[cymin] * double_f.y() * olap25 * zbkge / zgrid.width(czmin);
        drag_s(ix,cymax,czmin).zz() += yoverlap[cymax] * double_f.y() * olap25 * zbkge / zgrid.width(czmin);
        drag_s(ix,cymin,czmax).zz() += yoverlap[cymin] * double_f.y() * olap25 * zbkge / zgrid.width(czmax);
        drag_s(ix,cymax,czmax).zz() += yoverlap[cymax] * double_f.y() * olap25 * zbkge / zgrid.width(czmax);

        // Porous obstacles do not only have drag at edges
        if (xbkge < 1.0)
        {
            for (label iy = cymin+1; iy < cymax; iy++)
            {
                for (label iz = czmin+1; iz < czmax; iz++)
                {
                    // Internal only
                    drag_s(ix,iy,iz).xx() = xbkge / xgrid.width(ix);
                }
            }
        }
    }

    for (label iy = cymin; iy <= cymax; iy++)
    {
        scalar olap25 = 0.25 * yoverlap[iy];
        const scalar temp =
            ((pars.noIntersectN && vbkge < 0.334) ? 0 : (olap25 * vbkge));

        obs_count(cxmin,iy,czmin) += temp;
        obs_count(cxmax,iy,czmin) += temp;
        obs_count(cxmin,iy,czmax) += temp;
        obs_count(cxmax,iy,czmax) += temp;

        sub_count(cxmin,iy,czmin).y() += temp;
        sub_count(cxmax,iy,czmin).y() += temp;
        sub_count(cxmin,iy,czmax).y() += temp;
        sub_count(cxmax,iy,czmax).y() += temp;

        olap25 *= 2.0;

        if (iy > cymin && iy < cymax) // Avoid re-doing corners already done above
        {
            drag_s(cxmin,iy,czmin).zz() += xoverlap[cxmin] * double_f.x() * olap25 * zbkge / zgrid.width(czmin);
            drag_s(cxmax,iy,czmin).zz() += xoverlap[cxmin] * double_f.x() * olap25 * zbkge / zgrid.width(czmin);
            drag_s(cxmin,iy,czmax).zz() += xoverlap[cxmax] * double_f.x() * olap25 * zbkge / zgrid.width(czmax);
            drag_s(cxmax,iy,czmax).zz() += xoverlap[cxmax] * double_f.x() * olap25 * zbkge / zgrid.width(czmax);
        }
        drag_s(cxmin,iy,czmin).xx() += zoverlap[czmin] * double_f.z() * olap25 * xbkge / xgrid.width(cxmin);
        drag_s(cxmax,iy,czmin).xx() += zoverlap[czmax] * double_f.z() * olap25 * xbkge / xgrid.width(cxmax);
        drag_s(cxmin,iy,czmax).xx() += zoverlap[czmin] * double_f.z() * olap25 * xbkge / xgrid.width(cxmin);
        drag_s(cxmax,iy,czmax).xx() += zoverlap[czmax] * double_f.z() * olap25 * xbkge / xgrid.width(cxmax);

        // Porous obstacles do not only have drag at edges
        if (ybkge < 1.0)
        {
            for (label iz = czmin+1; iz < czmax; iz++)
            {
                for (label ix = cxmin+1; ix < cxmax; ix++)
                {
                    // Internal only
                    drag_s(ix,iy,iz).yy() = ybkge / ygrid.width(iy);
                }
            }
        }
    }

    for (label iz = czmin; iz <= czmax; iz++)
    {
        scalar olap25 = 0.25 * zoverlap[iz];
        const scalar temp =
            ((pars.noIntersectN && vbkge < 0.334) ? 0 : (olap25 * vbkge));

        obs_count(cxmin,cymin,iz) += temp;
        obs_count(cxmin,cymax,iz) += temp;
        obs_count(cxmax,cymin,iz) += temp;
        obs_count(cxmax,cymax,iz) += temp;

        sub_count(cxmin,cymin,iz).z() += temp;
        sub_count(cxmin,cymax,iz).z() += temp;
        sub_count(cxmax,cymin,iz).z() += temp;
        sub_count(cxmax,cymax,iz).z() += temp;

        olap25 *= 2.0;

        if (iz > czmin && iz < czmax) // Avoid re-doing corners already done above
        {
            drag_s(cxmin,cymin,iz).xx() += yoverlap[cymin] * double_f.y() * olap25 * xbkge / xgrid.width(cxmin);
            drag_s(cxmax,cymin,iz).xx() += yoverlap[cymin] * double_f.y() * olap25 * xbkge / xgrid.width(cxmax);
            drag_s(cxmin,cymax,iz).xx() += yoverlap[cymax] * double_f.y() * olap25 * xbkge / xgrid.width(cxmin);
            drag_s(cxmax,cymax,iz).xx() += yoverlap[cymax] * double_f.y() * olap25 * xbkge / xgrid.width(cxmax);

            drag_s(cxmin,cymin,iz).yy() += xoverlap[cxmin] * double_f.x() * olap25 * ybkge / ygrid.width(cymin);
            drag_s(cxmax,cymin,iz).yy() += xoverlap[cxmax] * double_f.x() * olap25 * ybkge / ygrid.width(cymin);
            drag_s(cxmin,cymax,iz).yy() += xoverlap[cxmin] * double_f.x() * olap25 * ybkge / ygrid.width(cymax);
            drag_s(cxmax,cymax,iz).yy() += xoverlap[cxmax] * double_f.x() * olap25 * ybkge / ygrid.width(cymax);
        }

        // Porous obstacles do not only have drag at edges
        if (zbkge < 1.0)
        {
            for (label ix = cxmin+1; ix < cxmax; ix++)
            {
                for (label iy = cymin+1; iy < cymax; iy++)
                {
                    // Internal only
                    drag_s(ix,iy,iz).zz() = zbkge / zgrid.width(iz);
                }
            }
        }
    }
}


// ************************************************************************* //
