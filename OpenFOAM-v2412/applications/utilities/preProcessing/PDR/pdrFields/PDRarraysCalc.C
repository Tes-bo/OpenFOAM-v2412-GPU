/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "PDRarrays.H"
#include "PDRblock.H"
#include "PDRpatchDef.H"
#include "PDRmeshArrays.H"
#include "PDRparams.H"

#include "PDRsetFields.H"

#include "bitSet.H"
#include "DynamicList.H"
#include "dimensionSet.H"
#include "symmTensor.H"
#include "SquareMatrix.H"
#include "IjkField.H"
#include "MinMax.H"
#include "volFields.H"
#include "OFstream.H"
#include "OSspecific.H"

#ifndef FULLDEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// A good ijk index has all components >= 0
static inline bool isGoodIndex(const Foam::labelVector& idx)
{
    return (idx.x() >= 0 && idx.y() >= 0 && idx.z() >= 0);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

static Foam::HashTable<Foam::string> fieldNotes
({
    { "Lobs", "" },
    { "Aw", "surface area per unit volume" },
    { "CR", "" },
    { "CT", "" },
    { "N", "" },
    { "ns", "" },
    { "Nv", "" },
    { "nsv", "" },
    { "Bv", "area blockage" },
    { "B", "area blockage" },
    { "betai", "" },
    { "Blong", "longitudinal blockage" },
    { "Ep", "1/Lobs" },
});


// calc_fields


// Local Functions
/*
// calc_drag_etc
make_header
tail_field
write_scalarField
write_uniformField
write_symmTensorField
write_pU_fields
write_blocked_face_list
write_blockedCellsSet
*/

// Somewhat similar to what the C-fprintf would have had
static constexpr unsigned outputPrecision = 8;

void calc_drag_etc
(
    double  brs, double  brr, bool blocked,
    double  surr_br, double  surr_dr,
    scalar* drags_p, scalar* dragr_p,
    double  count,
    scalar* cbdi_p,
    double cell_vol
);


void write_scalarField
(
    const word& fieldName, const IjkField<scalar>& fld,
    const scalar& deflt, const scalarMinMax& limits, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
);

void write_uniformField
(
    const word& fieldName, const scalar& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
);

void write_pU_fields
(
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const fileName& casepath
);

void write_symmTensorField
(
    const word& fieldName, const IjkField<symmTensor>& fld,
    const symmTensor& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
);

void write_symmTensorFieldV
(
    const word& fieldName, const IjkField<vector>& fld,
    const symmTensor& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
);

void write_blocked_face_list
(
    const IjkField<vector>& face_block,
    const IjkField<labelVector>& face_patch,
    const IjkField<scalar>& obs_count,
    IjkField<vector>& sub_count,
    IjkField<Vector<direction>>& n_blocked_faces,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    double limit_par, const fileName& casepath
);

void write_blockedCellsSet
(
    const IjkField<scalar>& fld,
    const PDRmeshArrays& meshIndexing, double limit_par,
    const IjkField<Vector<direction>>& n_blocked_faces,
    const fileName& casepath,
    const word& listName
);


// The average values of surrounding an array position
static inline scalar averageSurrounding
(
    const SquareMatrix<scalar>& mat,
    const label i,
    const label j
)
{
    return
    (
        mat(i,j)   + mat(i,j+1) + mat(i,j+2)
      + mat(i+1,j) /* centre */ + mat(i+1,j+2)
      + mat(i+2,j) + mat(i+2,j+1) + mat(i+2,j+2)
    ) / 8.0;  // Average
}


// Helper
template<class Type>
static inline Ostream& putUniform(Ostream& os, const word& key, const Type& val)
{
    os.writeKeyword(key);
    os  << word("uniform") << token::SPACE << val;
    os.endEntry();
    return os;
}


static void make_header
(
    Ostream& os,
    const fileName& location,
    const word& clsName,
    const word& object
)
{
    string note = fieldNotes(object);

    IOobject::writeBanner(os);

    os  << "FoamFile\n{\n"
        << "    version     2.0;\n"
        << "    format      ascii;\n"
        << "    class       " << clsName << ";\n";

    if (!note.empty())
    {
        os << "    note        " << note << ";\n";
    }

    if (!location.empty())
    {
        os << "    location    " << location << ";\n";
    }

    os  << "    object      " << object << ";\n"
        << "}\n";

    IOobject::writeDivider(os) << nl;
}


void Foam::PDRarrays::calculateAndWrite
(
    PDRarrays& arr,
    const PDRmeshArrays& meshIndexing,
    const fileName& casepath,
    const UList<PDRpatchDef>& patches
)
{
    if (isNull(arr.block()))
    {
        FatalErrorInFunction
            << "No PDRblock set" << nl
            << exit(FatalError);
    }

    const PDRblock& pdrBlock = arr.block();

    const labelVector& cellDims = meshIndexing.cellDims;
    const labelVector& faceDims = meshIndexing.faceDims;

    const int xdim = faceDims.x();
    const int ydim = faceDims.y();
    const int zdim = faceDims.z();
    const scalar maxCT = pars.maxCR * pars.cb_r;


    // Later used to store the total effective blockage ratio per cell/direction
    IjkField<symmTensor>& drag_s = arr.drag_s;

    IjkField<vector>& drag_r = arr.drag_r;

    const IjkField<vector>& area_block_s = arr.area_block_s;
    const IjkField<vector>& area_block_r = arr.area_block_r;
    const IjkField<Vector<bool>>& dirn_block = arr.dirn_block;

    const IjkField<vector>& betai_inv1 = arr.betai_inv1;

    IjkField<scalar>& obs_count = arr.obs_count;
    IjkField<vector>& sub_count = arr.sub_count; // ns. Later used to hold longitudinal blockage
    const IjkField<vector>& grating_count = arr.grating_count;

    IjkField<scalar>& v_block = arr.v_block;
    IjkField<scalar>& surf = arr.surf;

    // Lobs. Later used for initial Ep
    IjkField<scalar>& obs_size = arr.obs_size;

    Info<< "Calculating fields" << nl;

    // Local scratch arrays

    // The turbulance generation field CT.
    // Later used to to hold the beta_i in tensor form
    IjkField<vector> cbdi(cellDims, Zero);


    // For 2D addressing it is convenient to just use the max dimension
    // and avoid resizing when handling each direction.

    // Dimension of the cells and a layer of surrounding halo cells
    const labelVector surrDims = (faceDims + labelVector::uniform(2));

    // Max addressing dimensions
    const label maxDim = cmptMax(surrDims);

    // Blockage-ratio correction to the drag
    //
    // neiBlock:
    //   2-D for averaging the blockage ratio of neighbouring cells.
    //   It extends one cell  outside the domain in each direction,
    //   so the indices are offset by 1.
    // neiDrag:
    //   2-D array for averaging the drag ratio of neighbouring cells

    SquareMatrix<scalar> neiBlock(maxDim, Zero);
    SquareMatrix<scalar> neiDrag(maxDim, Zero);

    // X blockage, drag

    for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
    {
        for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
        {
            for (label iz = 0; iz <= zdim; ++iz)
            {
                const label izz =
                    (iz == 0 ? 0 : iz == zdim ? zdim - 2 : iz - 1);

                neiBlock(iy+1, iz) =
                (
                    area_block_s(ix,iy,izz).x()
                  + area_block_r(ix,iy,izz).x()
                );

                neiDrag(iy+1, iz) =
                (
                    drag_s(ix,iy,izz).xx() * pars.cd_s
                  + drag_r(ix,iy,izz).x() * pars.cd_r
                );
            }
        }
        for (label iz = 0; iz < surrDims.z(); ++iz)
        {
            if (pars.yCyclic)
            {
                // Cyclic in y
                neiBlock(0, iz) = neiBlock(cellDims.y(), iz);
                neiDrag(0, iz)  = neiDrag(cellDims.y(), iz);
                neiBlock(ydim, iz) = neiBlock(1, iz);
                neiDrag(ydim, iz)  = neiDrag(1, iz);
            }
            else
            {
                neiBlock(0, iz) = neiBlock(1, iz);
                neiDrag(0, iz)  = neiDrag(1, iz);
                neiBlock(ydim, iz) = neiBlock(cellDims.y(), iz);
                neiDrag(ydim, iz)  = neiDrag(cellDims.y(), iz);
            }
        }

        for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
        {
            for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
            {
                const scalar cell_vol = pdrBlock.V(ix,iy,iz);

                const scalar surr_br = averageSurrounding(neiBlock, iy, iz);
                const scalar surr_dr = averageSurrounding(neiDrag, iy, iz);

                calc_drag_etc
                (
                    area_block_s(ix,iy,iz).x(),
                    area_block_r(ix,iy,iz).x(),
                    dirn_block(ix,iy,iz).x(),
                    surr_br, surr_dr,
                    &(drag_s(ix,iy,iz).xx()),
                    &(drag_r(ix,iy,iz).x()),
                    obs_count(ix,iy,iz),
                    &(cbdi(ix,iy,iz).x()),
                    cell_vol
                );
            }
        }
    }


    // Y blockage, drag

    neiBlock = Zero;
    neiDrag = Zero;

    for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
    {
        for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
        {
            for (label ix = 0; ix <= xdim; ++ix)
            {
                const label ixx =
                    (ix == 0 ? 0 : ix == xdim ? xdim - 2 : ix - 1);

                neiBlock(iz+1, ix) =
                (
                    area_block_s(ixx,iy,iz).y()
                  + area_block_r(ixx,iy,iz).y()
                );
                neiDrag(iz+1, ix) =
                (
                    drag_s(ixx,iy,iz).yy() * pars.cd_s
                  + drag_r(ixx,iy,iz).y() * pars.cd_r
                );
            }
        }
        for (label ix = 0; ix < surrDims.x(); ++ix)
        {
            neiBlock(0, ix) = neiBlock(1, ix);
            neiDrag(0, ix)  = neiDrag(1, ix);
            neiBlock(zdim, ix) = neiBlock(cellDims.z(), ix);
            neiDrag(zdim, ix)  = neiDrag(cellDims.z(), ix);
        }

        for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
        {
            for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
            {
                const scalar cell_vol = pdrBlock.V(ix,iy,iz);

                const scalar surr_br = averageSurrounding(neiBlock, iz, ix);
                const scalar surr_dr = averageSurrounding(neiDrag, iz, ix);

                calc_drag_etc
                (
                    area_block_s(ix,iy,iz).y(),
                    area_block_r(ix,iy,iz).y(),
                    dirn_block(ix,iy,iz).y(),
                    surr_br, surr_dr,
                    &(drag_s(ix,iy,iz).yy()),
                    &(drag_r(ix,iy,iz).y()),
                    obs_count(ix,iy,iz),
                    &(cbdi(ix,iy,iz).y()),
                    cell_vol
                );
            }
        }
    }


    // Z blockage, drag

    neiBlock = Zero;
    neiDrag = Zero;

    for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
    {
        for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
        {
            for (label iy = 0; iy <= ydim; ++iy)
            {
                label iyy;

                if (pars.yCyclic)
                {
                    iyy = (iy == 0 ? ydim - 2 : iy == ydim ? 0 : iy - 1);
                }
                else
                {
                    iyy = (iy == 0 ? 0 : iy == ydim ? ydim - 2 : iy - 1);
                }

                neiBlock(ix+1, iy) =
                (
                    area_block_s(ix,iyy,iz).z()
                  + area_block_r(ix,iyy,iz).z()
                );
                neiDrag(ix+1, iy) =
                (
                    drag_s(ix,iyy,iz).zz() * pars.cd_s
                  + drag_r(ix,iyy,iz).z() * pars.cd_r
                );
            }
        }
        for (label iy = 0; iy < surrDims.y(); ++iy)
        {
            neiBlock(0, iy) = neiBlock(1, iy);
            neiDrag(0, iy)  = neiDrag(1, iy);
            neiBlock(xdim, iy) = neiBlock(cellDims.x(), iy);
            neiDrag(xdim, iy)  = neiDrag(cellDims.x(), iy);
        }

        for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
        {
            for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
            {
                const scalar cell_vol = pdrBlock.V(ix,iy,iz);

                const scalar surr_br = averageSurrounding(neiBlock, ix, iy);
                const scalar surr_dr = averageSurrounding(neiDrag, ix, iy);

                calc_drag_etc
                (
                    area_block_s(ix,iy,iz).z(),
                    area_block_r(ix,iy,iz).z(),
                    dirn_block(ix,iy,iz).z(),
                    surr_br, surr_dr,
                    &(drag_s(ix,iy,iz).zz()),
                    &(drag_r(ix,iy,iz).z()),
                    obs_count(ix,iy,iz),
                    &(cbdi(ix,iy,iz).z()),
                    cell_vol
                );
            }
        }
    }

    neiBlock.clear();
    neiDrag.clear();


    // Calculate other parameters

    for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
    {
        for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
        {
            for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
            {
                const scalar dx = pdrBlock.dx(ix);
                const scalar dy = pdrBlock.dy(iy);
                const scalar dz = pdrBlock.dz(iz);
                const scalar cell_vol = pdrBlock.V(ix, iy, iz);
                const scalar cell_size = pdrBlock.width(ix, iy, iz);

                drag_s(ix,iy,iz).xy() *= pars.cd_s;
                drag_s(ix,iy,iz).xz() *= pars.cd_s;
                drag_s(ix,iy,iz).yz() *= pars.cd_s;

                if (drag_s(ix,iy,iz).xx() > pars.maxCR) { drag_s(ix,iy,iz).xx() = pars.maxCR; } ;
                if (drag_s(ix,iy,iz).yy() > pars.maxCR) { drag_s(ix,iy,iz).yy() = pars.maxCR; } ;
                if (drag_s(ix,iy,iz).zz() > pars.maxCR) { drag_s(ix,iy,iz).zz() = pars.maxCR; } ;

                if (cbdi(ix,iy,iz).x() > maxCT ) { cbdi(ix,iy,iz).x() = maxCT; } ;
                if (cbdi(ix,iy,iz).y() > maxCT ) { cbdi(ix,iy,iz).y() = maxCT; } ;
                if (cbdi(ix,iy,iz).z() > maxCT ) { cbdi(ix,iy,iz).z() = maxCT; } ;

                surf(ix,iy,iz) /= cell_vol;

                /* Calculate length scale of obstacles in each cell
                 Result is stored in surf. */

                {
                    const scalar vb = v_block(ix,iy,iz);

                    if
                    (
                        (
                            ((area_block_s(ix,iy,iz).x() + area_block_r(ix,iy,iz).x()) < MIN_AB_FOR_SIZE)
                         && ((area_block_s(ix,iy,iz).y() + area_block_r(ix,iy,iz).y()) < MIN_AB_FOR_SIZE)
                         && ((area_block_s(ix,iy,iz).z() + area_block_r(ix,iy,iz).z()) < MIN_AB_FOR_SIZE)
                        )
                     || ( vb > MAX_VB_FOR_SIZE )
                     || ((obs_count(ix,iy,iz) + cmptSum(grating_count(ix,iy,iz))) < MIN_COUNT_FOR_SIZE)
                     || ( surf(ix,iy,iz) <= 0.0 )
                    )
                    {
                        obs_size(ix,iy,iz) = cell_size * pars.empty_lobs_fac;
                    }
                    else
                    {
                        /* A small sliver of a large cylinder ina cell can give large surface area
                         but low volume, hence snall "size". Therefore the vol/area formulation
                         is only fully implemented when count is at least COUNT_FOR_SIZE.*/
                        double nn, lobs, lobsMax;
                        nn = obs_count(ix,iy,iz) - sub_count(ix,iy,iz).x() + grating_count(ix,iy,iz).x();
                        if ( nn < 1.0 ) { nn = 1.0; }
                        lobsMax = (area_block_s(ix,iy,iz).x() + area_block_r(ix,iy,iz).x()) / nn * std::sqrt( dy * dz );
                        nn = obs_count(ix,iy,iz) - sub_count(ix,iy,iz).y() + grating_count(ix,iy,iz).y();
                        if ( nn < 1.0 ) { nn = 1.0; }
                        lobs = (area_block_s(ix,iy,iz).y() + area_block_r(ix,iy,iz).y()) / nn * std::sqrt( dz * dx );
                        if ( lobs > lobsMax )
                        {
                            lobsMax = lobs;
                        }

                        nn = obs_count(ix,iy,iz) - sub_count(ix,iy,iz).z() + grating_count(ix,iy,iz).z();
                        if ( nn < 1.0 ) { nn = 1.0; }
                        lobs = (area_block_s(ix,iy,iz).z() + area_block_r(ix,iy,iz).z()) / nn * std::sqrt( dx * dy );
                        if ( lobs > lobsMax )
                        {
                            lobsMax = lobs;
                        }

                        obs_size(ix,iy,iz) = lobsMax;
                    }
                }

                /* The formulation correctly deals with triple intersections. For quadruple intersections
                 and worse, there are very many second level overlaps and the resulting volume can be large
                 positive. However, many or all of these may be eliminated because of the minimum volume of
                 overlap blocks. Then the result can be negative volume - constrain accordingly
                 */

                if (v_block(ix,iy,iz) < 0)
                {
                    v_block(ix,iy,iz) = 0;
                }
                else if (v_block(ix,iy,iz) > 1)
                {
                    v_block(ix,iy,iz) = 1;
                }

                /* We can get -ve sub_count (ns) if two pipes/bars intersect and the dominat direction
                 of the (-ve) intersection block is not the same as either of the intersecting obstacles.
                 Also, if we have two hirizontal abrs intersecting, the overlap block can have vertical
                 edges in a cell where the original bars do not. This can give -ve N and ns.
                 Negative N is removed by  write_scalar.  */

                for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
                {
                    if (sub_count(ix,iy,iz)[cmpt] < 0)
                    {
                        sub_count(ix,iy,iz)[cmpt] = 0;
                    }
                }

                v_block(ix,iy,iz) = 1.0 - v_block(ix,iy,iz); // Now porosity
            }
        }
    }


//*** Now we start writing the fields *********//

        /* v_block is now porosity
           The maximum value does not override the default value placed in the external cells,
           so pars.cong_max_betav can be set just below 1 to mark the congested-region cells
           for use by the adaptive mesh refinement.   */

    IjkField<Vector<direction>> n_blocked_faces
    (
        faceDims,
        Vector<direction>::uniform(0)
    );

    write_blocked_face_list
    (
        arr.face_block, arr.face_patch,
        obs_count, sub_count, n_blocked_faces,
        meshIndexing, patches,
        pars.blockedFacePar, casepath
    );
    write_blockedCellsSet
    (
        arr.v_block,
        meshIndexing, pars.blockedCellPoros, n_blocked_faces,
        casepath, "blockedCellsSet"
    );

    write_scalarField
    (
        "betav", arr.v_block, 1, {0, pars.cong_max_betav}, "zeroGradient",
        meshIndexing, patches,
        dimless, casepath
    );

    for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
    {
        for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
        {
            for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
            {
                const scalar cell_vol = pdrBlock.V(ix, iy, iz);

                /* After the correction to set the number of obstacles normal to a blocked face
                 to be zero, we can have N and all the components of ns the same. Then there
                 are no obstacles in the cell as the number in each direction is n minus ns component),
                 but N is not zero. This can cause problems. We reduce all four numbers by the same amount,
                 which is OK as only the difference is used except when N is checked to se if there are
                 any obstacles in then cell. */

                scalar nmin = cmptMin(sub_count(ix,iy,iz));

                sub_count(ix,iy,iz).x() -= nmin;
                sub_count(ix,iy,iz).y() -= nmin;
                sub_count(ix,iy,iz).z() -= nmin;

                obs_count(ix,iy,iz) -= nmin;

                assert(obs_count(ix,iy,iz) > -1);
                if ( pars.new_fields )
                {
                    /* New fields Nv and nsv are intensive quantities that stay unchanged as a cell is subdivided
                     We do not divide by cell volume because we assume that typical obstacle
                     is a cylinder passing through the cell */
                    const scalar cell_23 = ::pow(cell_vol, 2.0/3.0);
                    obs_count(ix,iy,iz) /= cell_23;
                    sub_count(ix,iy,iz) /= cell_23;
                }
            }
        }
    }


    {
        Info<< "Writing field files" << nl;

        // obs_size is now the integral scale of the generated turbulence
        write_scalarField
        (
            "Lobs", arr.obs_size, DEFAULT_LOBS, {0, 10}, "zeroGradient",
            meshIndexing, patches,
            dimLength, casepath
        );
        // surf is now surface area per unit volume
        write_scalarField
        (
            "Aw", arr.surf, 0, {0, 1000}, "zeroGradient",
            meshIndexing, patches,
            inv(dimLength), casepath
        );
        write_symmTensorField
        (
            "CR", arr.drag_s, Zero, "zeroGradient",
            meshIndexing, patches, inv(dimLength), casepath
        );
        write_symmTensorFieldV
        (
            "CT", cbdi, Zero, "zeroGradient",
            meshIndexing, patches,
            inv(dimLength), casepath
        );
        if ( pars.new_fields )
        {
            // These have been divided by cell volume ^ (2/3)
            write_scalarField
            (
                "Nv", arr.obs_count, 0, {0, 1000}, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_symmTensorFieldV
            (
                "nsv", arr.sub_count, Zero, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }
        else
        {
            write_scalarField
            (
                "N", arr.obs_count, 0, {0, 1000}, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_symmTensorFieldV
            (
                "ns", arr.sub_count, Zero, "zeroGradient",
                meshIndexing, patches, dimless, casepath
            );
        }

        // Compute some further variables; store in already used arrays
        // Re-use the drag array
        drag_s = Zero;

        for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
        {
            for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
            {
                for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
                {
                    // Effective blockage ratio per cell/direction
                    vector eff_block =
                    (
                        area_block_s(ix,iy,iz) * pars.cd_s/pars.cd_r
                      + area_block_r(ix,iy,iz)
                    );

                    // Convert from B to Bv
                    if (pars.new_fields)
                    {
                        eff_block /= pdrBlock.width(ix, iy, iz);
                    }

                    // Effective blockage is zero when faces are blocked
                    for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
                    {
                        if (dirn_block(ix,iy,iz)[cmpt] || eff_block[cmpt] < 0)
                        {
                            eff_block[cmpt] = 0;
                        }
                    }

                    // Use the drag array to store the total effective blockage ratio per cell/direction
                    // - off-diagonal already zeroed
                    drag_s(ix,iy,iz).xx() = eff_block.x();
                    drag_s(ix,iy,iz).yy() = eff_block.y();
                    drag_s(ix,iy,iz).zz() = eff_block.z();

                    cbdi(ix,iy,iz).x() = 1.0 / (betai_inv1(ix,iy,iz).x() + 1.0);
                    cbdi(ix,iy,iz).y() = 1.0 / (betai_inv1(ix,iy,iz).y() + 1.0);
                    cbdi(ix,iy,iz).z() = 1.0 / (betai_inv1(ix,iy,iz).z() + 1.0);

                    if (cbdi(ix,iy,iz).z() < 0 || cbdi(ix,iy,iz).z() > 1.0)
                    {
                        WarningInFunction
                            << "beta_i problem. z-betai_inv1=" << betai_inv1(ix,iy,iz).z()
                            << " beta_i=" << cbdi(ix,iy,iz).z()
                            << nl;
                    }

                    //Use the obs_size array to store Ep
                    //We use Ep/(Xp-0.999) as length scale to avoid divide by zero,
                    // so this is OK for initial Xp=1.
                    obs_size(ix,iy,iz) = 0.001 / obs_size(ix,iy,iz);

                    // Use the count array to store the combustion flag ( --1 everywhere in rectangular cells).
                    obs_count(ix,iy,iz) = 1.0;
                }
            }
        }

        // drag array holds area blockage
        if ( pars.new_fields )
        {
            write_symmTensorField
            (
                "Bv", arr.drag_s, Zero, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }
        else
        {
            write_symmTensorField
            (
                "B", arr.drag_s, Zero, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }

        // cbdi array holds beta_i
        write_symmTensorFieldV
        (
            "betai", cbdi, symmTensor::I, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );

        // The longitudinal blockage
        write_symmTensorFieldV
        (
            "Blong", arr.along_block, Zero, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );

        // obs_size array now contains 1/Lobs
        write_scalarField
        (
            "Ep", arr.obs_size, DEFAULT_EP, {0, 10}, "zeroGradient",
            meshIndexing, patches,
            inv(dimLength), casepath
        );
        write_uniformField
        (
            "b", 1.0, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );
        write_uniformField
        (
            "k", DEFAULT_K, K_WALL_FN,
            meshIndexing, patches,
            sqr(dimVelocity),
            casepath
        );

        write_uniformField
        (
            "epsilon", DEFAULT_EPS, EPS_WALL_FN,
            meshIndexing, patches,
            sqr(dimVelocity)/dimTime, casepath
        );
        write_uniformField
        (
            "ft", 0, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );
        write_uniformField
        (
            "Su", DEFAULT_SU, "zeroGradient",
            meshIndexing, patches,
            dimVelocity, casepath
        );
        write_uniformField
        (
            "T", DEFAULT_T, "zeroGradient",
            meshIndexing, patches,
            dimTemperature, casepath
        );
        write_uniformField
        (
            "Tu", DEFAULT_T, "zeroGradient",
            meshIndexing,  patches,
            dimTemperature, casepath
        );
        write_uniformField
        (
            "Xi", 1, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );
        write_uniformField
        (
            "Xp", 1, "zeroGradient",
            meshIndexing, patches,
            dimless, casepath
        );
        write_uniformField
        (
            "GRxp", 0, "zeroGradient",
            meshIndexing, patches,
            inv(dimTime), casepath
        );
        write_uniformField
        (
            "GRep", 0, "zeroGradient",
            meshIndexing, patches,
            inv(dimLength*dimTime), casepath
        );
        write_uniformField
        (
            "RPers", 0, "zeroGradient",
            meshIndexing, patches,
            inv(dimTime), casepath
        );
        write_pU_fields(meshIndexing, patches, casepath);

        write_uniformField
        (
            "alphat", 0, ALPHAT_WALL,
            meshIndexing, patches,
            dimMass/(dimLength*dimTime),
            casepath
        );
        write_uniformField
        (
            "nut", 0, NUT_WALL_FN,
            meshIndexing, patches,
            dimViscosity, casepath
        );
        // combustFlag is 1 in rectangular region, 0 or 1 elsewhere
        // (although user could set it to another value)
        if (equal(pars.outerCombFac, 1))
        {
            write_uniformField
            (
                "combustFlag", 1, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }
        else
        {
            write_scalarField
            (
                "combustFlag", arr.obs_count, pars.outerCombFac, {0, 1}, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }
        if ( pars.deluge )
        {
            write_uniformField
            (
                "H2OPS", 0, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_uniformField
            (
                "AIR", 0, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_uniformField
            (
                "Ydefault", 0, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_uniformField
            (
                "eRatio", 1, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
            write_uniformField
            (
                "sprayFlag", 1, "zeroGradient",
                meshIndexing, patches,
                dimless, casepath
            );
        }
    }
}


void Foam::PDRarrays::calculateAndWrite
(
    const fileName& casepath,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches
)
{
    calculateAndWrite(*this, meshIndexing, casepath, patches);
}


void calc_drag_etc
(
    double brs, double  brr, bool  blocked,
    double surr_br, double surr_dr,
    scalar* drags_p, scalar* dragr_p,
    double count,
    scalar* cbdi_p,
    double cell_vol
)
{
    // Total blockage ratio
    scalar br = brr + brs;

    // Idealise obstacle arrangement as sqrt(count) rows.
    // Make br the blockage ratio for each row.
    if (count > 1.0) { br /= std::sqrt(count); }

    const scalar alpha =
    (
        br < 0.99
      ? (1.0 - 0.5 * br) / (1.0 - br) / (1.0 - br)
      : GREAT
    );

    // For the moment keep separate the two contributions to the blockage-corrected drag
    /* An isolated long obstcale will have two of the surronding eight cells with the same blockage,
     so surr_br would be br/4. In this case no correction. Rising to full correction when
     all surrounding cells have the same blockage. */
    const scalar expon =
    (
        br > 0.0
      ? min(max((surr_br / br - 0.25) * 4.0 / 3.0, scalar(0)), scalar(1))
      : 0.0
    );

    const scalar alpha_r = ::pow(alpha, 0.5 + 0.5 * expon);
    const scalar alpha_s = ::pow(alpha, expon);

    *dragr_p *=  alpha_r;
    *drags_p *=  ::pow(alpha_s, 1.09);
    *cbdi_p  = ( pars.cb_r * pars.cd_r * *dragr_p + pars.cb_s * pars.cd_s * *drags_p );
    if ( *cbdi_p < 0.0 ) { *cbdi_p = 0.0; }

    // Finally sum the drag.
    *drags_p = ( *drags_p * pars.cd_s + *dragr_p * pars.cd_r );
    if ( *drags_p < 0.0 ) { *drags_p = 0.0; }
    /* If well-blocked cells are surrounded by empty cells, the flow just goes round
     and the drag parameters have little effect. So, for any cells much more empty
     than the surrounding cells, we put some CR in there as well.                  */
    if ( (surr_dr * 0.25) > *drags_p )
    {
        *drags_p = surr_dr * 0.25;
        *cbdi_p = *drags_p * (pars.cb_r + pars.cb_s ) * 0.5;
        // Don't know whether surr. stuff was round or sharp; use average of cb factors
    }
    if ( blocked ) { *cbdi_p = 0.0; *drags_p = 0.0; *dragr_p = 0.0; }
}


void Foam::PDRarrays::blockageSummary() const
{
    if (isNull(block()))
    {
        WarningInFunction
            << nl
            << "No blockage information - PDRblock is not set" << nl;
        return;
    }

    const PDRblock& pdrBlock = block();

    scalar totArea = 0;
    scalar totCount = 0;
    scalar totVolBlock = 0;

    vector totBlock(Zero);
    vector totDrag(Zero);

    for (label iz = 0; iz < pdrBlock.size(vector::Z); ++iz)
    {
        for (label iy = 0; iy < pdrBlock.size(vector::Y); ++iy)
        {
            for (label ix = 0; ix < pdrBlock.size(vector::X); ++ix)
            {
                const labelVector ijk(ix,iy,iz);

                totVolBlock += v_block(ijk) * pdrBlock.V(ijk);
                totArea += surf(ijk);

                totCount += max(0, obs_count(ijk));

                totDrag.x() += max(0, drag_s(ijk).xx());
                totDrag.y() += max(0, drag_s(ijk).yy());
                totDrag.z() += max(0, drag_s(ijk).zz());

                for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
                {
                    totBlock[cmpt] += max(0, area_block_s(ijk)[cmpt]);
                    totBlock[cmpt] += max(0, area_block_r(ijk)[cmpt]);
                }
            }
        }
    }

    Info<< nl
        << "Volume blockage: " << totVolBlock << nl
        << "Total drag:  " << totDrag << nl
        << "Total count: " << totCount << nl
        << "Total area blockage: " << totBlock << nl
        << "Total surface area: " << totArea << nl;
}


// ------------------------------------------------------------------------- //

// Another temporary measure
template<class Type>
static void tail_field
(
    Ostream& os,
    const Type& deflt,
    const char* wall_bc,
    const UList<PDRpatchDef>& patches
)
{
    // ground
    {
        os.beginBlock(pars.groundPatchName);
        os.writeEntry("type", wall_bc);
        putUniform(os, "value", deflt);
        os.endBlock();
    }

    forAll(patches, patchi)
    {
        const word& patchName = patches[patchi].patchName;

        if (PDRpatchDef::BLOCKED_FACE == patchi)
        {
            // blockedFaces
            os.beginBlock(patchName);

            // No wall functions for blockedFaces patch unless selected
            if (pars.blockedFacesWallFn)
            {
                os.writeEntry("type", wall_bc);
                putUniform(os, "value", deflt);
            }
            else
            {
                os.writeEntry("type", "zeroGradient");
            }

            os.endBlock();
        }
        else if (patches[patchi].patchType == 0)
        {
            os.beginBlock(patchName);

            os.writeEntry("type", wall_bc);
            putUniform(os, "value", deflt);

            os.endBlock();
        }
        else
        {
            os.beginBlock(word(patchName + "Wall"));
            os.writeEntry("type", wall_bc);
            putUniform(os, "value", deflt);
            os.endBlock();

            os.beginBlock(word(patchName + "Cyclic_half0"));
            os.writeEntry("type", "cyclic");
            os.endBlock();

            os.beginBlock(word(patchName + "Cyclic_half1"));
            os.writeEntry("type", "cyclic");
            os.endBlock();
        }
    }

    if (pars.yCyclic)
    {
        os.beginBlock("Cyclic_half0");
        os.writeEntry("type", "cyclic");
        os.endBlock();

        os.beginBlock("Cyclic_half1");
        os.writeEntry("type", "cyclic");
        os.endBlock();
    }
    else
    {
        os.beginBlock("ySymmetry");
        os.writeEntry("type", "symmetryPlane");
        os.endBlock();
    }

    if (pars.two_d)
    {
        os.beginBlock("z_boundaries");
        os.writeEntry("type", "empty");
        os.endBlock();
    }

    if (pars.outer_orthog)
    {
        os.beginBlock("outer_inner");
        os.writeEntry("type", "cyclicAMI");
        os.writeEntry("neighbourPatch", "inner_outer");
        os.endBlock();

        os.beginBlock("inner_outer");
        os.writeEntry("type", "cyclicAMI");
        os.writeEntry("neighbourPatch", "outer_inner");
        os.endBlock();
    }
}


// ------------------------------------------------------------------------- //

void write_scalarField
(
    const word& fieldName, const IjkField<scalar>& fld,
    const scalar& deflt, const scalarMinMax& limits, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
)
{
    fileName path = (casepath / pars.timeName / fieldName);
    OFstream os(path);
    os.precision(outputPrecision);

    make_header(os, "", volScalarField::typeName, fieldName);

    os.writeEntry("dimensions", dims);

    os << nl;
    os.writeKeyword("internalField")
        << word("nonuniform") << token::SPACE << word("List<scalar>") << nl
        << meshIndexing.nCells() << nl << token::BEGIN_LIST << nl;

    for (label celli=0; celli < meshIndexing.nCells(); ++celli)
    {
        const labelVector& cellIdx = meshIndexing.cellIndex[celli];

        if (!isGoodIndex(cellIdx))
        {
            os  << deflt << nl;
            continue;
        }

        os << limits.clamp(fld(cellIdx)) << nl;
    }

    os << token::END_LIST;
    os.endEntry();

    os << nl;
    os.beginBlock("boundaryField");


    // outer
    {
        os.beginBlock(pars.outerPatchName);

        os.writeEntry("type", "inletOutlet");
        putUniform(os, "inletValue", deflt);
        putUniform(os, "value", deflt);

        os.endBlock();
    }

    tail_field(os, deflt, wall_bc, patches);

    os.endBlock(); // boundaryField

    IOobject::writeEndDivider(os);
}


// ------------------------------------------------------------------------- //

void write_uniformField
(
    const word& fieldName, const scalar& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
)
{
    OFstream os(casepath / pars.timeName / fieldName);
    os.precision(outputPrecision);

    make_header(os, "", volScalarField::typeName, fieldName);

    os.writeEntry("dimensions", dims);

    os << nl;
    putUniform(os, "internalField", deflt);

    os << nl;
    os.beginBlock("boundaryField");

    // outer
    {
        os.beginBlock(pars.outerPatchName);

        if (fieldName == "alphat" || fieldName == "nut")
        {
            // Different b.c. for alphat & nut
            os.writeEntry("type", "calculated");
        }
        else
        {
            os.writeEntry("type", "inletOutlet");
            putUniform(os, "inletValue", deflt);
        }

        putUniform(os, "value", deflt);
        os.endBlock();
    }

    tail_field(os, deflt, wall_bc, patches);

    os.endBlock(); // boundaryField

    IOobject::writeEndDivider(os);
}


// ------------------------------------------------------------------------- //

void write_pU_fields
(
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const fileName& casepath
)
{
    // Velocity field
    {
        OFstream os(casepath / pars.timeName / "U");
        os.precision(outputPrecision);

        make_header(os, "", volVectorField::typeName, "U");

        os.writeEntry("dimensions", dimVelocity);

        os << nl;
        putUniform(os, "internalField", vector::zero);

        os << nl;
        os.beginBlock("boundaryField");

        // outer
        {
            os.beginBlock(pars.outerPatchName);
            os.writeEntry("type", "inletOutlet");
            putUniform(os, "inletValue", vector::zero);
            os.endBlock();
        }

        // ground
        {
            os.beginBlock(pars.groundPatchName);
            os.writeEntry("type", "zeroGradient");
            os.endBlock();
        }

        // Patch 0 is the blocked faces' and 1 is mergingFaces for ignition cell
        for (label patchi = 0; patchi < 3; ++patchi)
        {
            os.beginBlock(patches[patchi].patchName);
            os.writeEntry("type", pars.UPatchBc.c_str());
            os.endBlock();
        }

        for (label patchi = 3; patchi < patches.size(); ++patchi)
        {
            const PDRpatchDef& p = patches[patchi];
            const word& patchName = p.patchName;

            if (p.patchType == 0)
            {
                os.beginBlock(patchName);

                os.writeEntry("type", "timeVaryingMappedFixedValue");
                os.writeEntry("fileName", "<case>" / (patchName + ".dat"));
                os.writeEntry("outOfBounds", "clamp");
                putUniform(os, "value", vector::zero);
                os.endBlock();
            }
            else
            {
                os.beginBlock(word(patchName + "Wall"));
                os.writeEntry("type", "activePressureForceBaffleVelocity");

                os.writeEntry("cyclicPatch", word(patchName + "Cyclic_half0"));
                os.writeEntry("openFraction", 0); // closed
                os.writeEntry("openingTime", p.blowoffTime);
                os.writeEntry("minThresholdValue", p.blowoffPress);
                os.writeEntry("maxOpenFractionDelta", 0.1);
                os.writeEntry("forceBased", "false");
                os.writeEntry("opening", "true");

                putUniform(os, "value", vector::zero);
                os.endBlock();

                os.beginBlock(word(patchName + "Cyclic_half0"));
                os.writeEntry("type", "cyclic");
                putUniform(os, "value", vector::zero);
                os.endBlock();

                os.beginBlock(word(patchName + "Cyclic_half1"));
                os.writeEntry("type", "cyclic");
                putUniform(os, "value", vector::zero);
                os.endBlock();
            }
        }

        if (pars.yCyclic)
        {
            os.beginBlock("yCyclic_half0");
            os.writeEntry("type", "cyclic");
            os.endBlock();

            os.beginBlock("yCyclic_half1");
            os.writeEntry("type", "cyclic");
            os.endBlock();
        }
        else
        {
            os.beginBlock("ySymmetry");
            os.writeEntry("type", "symmetryPlane");
            os.endBlock();
        }

        if ( pars.outer_orthog )
        {
            os.beginBlock("outer_inner");
            os.writeEntry("type", "cyclicAMI");
            os.writeEntry("neighbourPatch", "inner_outer");
            os.endBlock();

            os.beginBlock("inner_outer");
            os.writeEntry("type", "cyclicAMI");
            os.writeEntry("neighbourPatch", "outer_inner");
        }

        os.endBlock();  // boundaryField

        IOobject::writeEndDivider(os);
    }


    // Pressure field
    {
        const scalar deflt = DEFAULT_P;
        const char *wall_bc = "zeroGradient;\n\trho\trho";

        OFstream os(casepath / pars.timeName / "p");
        os.precision(outputPrecision);

        make_header(os, "", volScalarField::typeName, "p");

        os.writeEntry("dimensions", dimPressure);

        os << nl;
        putUniform(os, "internalField", deflt);

        os << nl;
        os.beginBlock("boundaryField");

        // outer
        {
            os.beginBlock(pars.outerPatchName);

            os.writeEntry("type", "waveTransmissive");
            os.writeEntry("gamma", 1.3);
            os.writeEntry("fieldInf", deflt);
            os.writeEntry("lInf", 5);
            putUniform(os, "value", deflt);
            os.endBlock();
        }

        tail_field(os, deflt, wall_bc, patches);

        os.endBlock();  // boundaryField

        IOobject::writeEndDivider(os);
    }
}


// ------------------------------------------------------------------------- //

void write_symmTensorField
(
    const word& fieldName,
    const IjkField<symmTensor>& fld,
    const symmTensor& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
)
{
    OFstream os(casepath / pars.timeName / fieldName);
    os.precision(outputPrecision);

    make_header(os, "", volSymmTensorField::typeName, fieldName);

    os.writeEntry("dimensions", dims);

    os << nl;
    os.writeKeyword("internalField")
        << word("nonuniform") << token::SPACE << word("List<symmTensor>") << nl
        << meshIndexing.nCells() << nl << token::BEGIN_LIST << nl;

    for (label celli=0; celli < meshIndexing.nCells(); ++celli)
    {
        const labelVector& cellIdx = meshIndexing.cellIndex[celli];

        if (!isGoodIndex(cellIdx))
        {
            os  << deflt << nl;
            continue;
        }

        os << fld(cellIdx) << nl;
    }
    os << token::END_LIST;
    os.endEntry();

    os << nl;
    os.beginBlock("boundaryField");

    // outer
    {
        os.beginBlock(pars.outerPatchName);

        os.writeEntry("type", "inletOutlet");
        putUniform(os, "inletValue", deflt);
        putUniform(os, "value", deflt);

        os.endBlock();
    }

    tail_field(os, deflt, wall_bc, patches);

    os.endBlock(); // boundaryField

    IOobject::writeEndDivider(os);
}


// Write a volSymmTensorField but with vectors as input.
// The off-diagonals are zero.
void write_symmTensorFieldV
(
    const word& fieldName,
    const IjkField<vector>& fld,
    const symmTensor& deflt, const char *wall_bc,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    const dimensionSet& dims, const fileName& casepath
)
{
    OFstream os(casepath / pars.timeName / fieldName);
    os.precision(outputPrecision);

    make_header(os, "", volSymmTensorField::typeName, fieldName);

    os.writeEntry("dimensions", dims);

    os << nl;
    os.writeKeyword("internalField")
        << word("nonuniform") << token::SPACE << word("List<symmTensor>") << nl
        << meshIndexing.nCells() << nl << token::BEGIN_LIST << nl;

    symmTensor val(symmTensor::zero);

    for (label celli=0; celli < meshIndexing.nCells(); ++celli)
    {
        const labelVector& cellIdx = meshIndexing.cellIndex[celli];

        if (!isGoodIndex(cellIdx))
        {
            os  << deflt << nl;
            continue;
        }

        const vector& vec = fld(cellIdx);

        val.xx() = vec.x();
        val.yy() = vec.y();
        val.zz() = vec.z();

        os << val << nl;
    }
    os << token::END_LIST;
    os.endEntry();

    os << nl;
    os.beginBlock("boundaryField");

    // outer
    {
        os.beginBlock(pars.outerPatchName);

        os.writeEntry("type", "inletOutlet");
        putUniform(os, "inletValue", deflt);
        putUniform(os, "value", deflt);

        os.endBlock();
    }

    tail_field(os, deflt, wall_bc, patches);

    os.endBlock(); // boundaryField

    IOobject::writeEndDivider(os);
}


// ------------------------------------------------------------------------- //

void write_blocked_face_list
(
    const IjkField<vector>& face_block,
    const IjkField<labelVector>& face_patch,
    const IjkField<scalar>& obs_count, IjkField<vector>& sub_count,
    IjkField<Vector<direction>>& n_blocked_faces,
    const PDRmeshArrays& meshIndexing,
    const UList<PDRpatchDef>& patches,
    double limit_par, const fileName& casepath
)
{
    /* Create the lists of face numbers for faces that have already been defined as
     belonging to (inlet) patches), and others that are found to be blocked.
     Then write these out to set files,    */

    const labelVector& cellDims = meshIndexing.cellDims;

    Map<bitSet> usedFaces;

    Info<< "Number of patches: " << patches.size() << nl;

    for (label facei=0; facei < meshIndexing.nFaces(); ++facei)
    {
        // The related i-j-k face index for the mesh face
        const labelVector& faceIdx = meshIndexing.faceIndex[facei];

        if (!isGoodIndex(faceIdx))
        {
            continue;
        }

        const label ix = faceIdx.x();
        const label iy = faceIdx.y();
        const label iz = faceIdx.z();
        const direction orient = meshIndexing.faceOrient[facei];

        label patchId = -1;
        scalar val(Zero);

        /* A bit messy to be changing sub_count here. but there is a problem of generation
         of subgrid flame area Xp when the flame approaches a blocked wall. the fix is to make
         the normal component of "n" zero in the cells adjacent to the blocked face. That component
         of n is zero when that component of sub_count i.e. ns) equals count (i.e. N). */
        {
            switch (orient)
            {
                case vector::X:
                {
                    // face_block is the face blockage;
                    // face_patch is the patch number on the face (if any)
                    val = face_block(faceIdx).x();
                    patchId = face_patch(faceIdx).x();

                    if
                    (
                        val > limit_par
                     && iy < cellDims[vector::Y]
                     && iz < cellDims[vector::Z]
                    )
                    {
                        // n_blocked_faces:
                        // count of x-faces blocked for this cell

                        if (ix < cellDims[vector::X])
                        {
                            ++n_blocked_faces(ix,iy,iz).x();
                            sub_count(ix,iy,iz).x() = obs_count(ix,iy,iz);
                        }

                        if (ix > 0)
                        {
                            // And the neighbouring cell
                            ++n_blocked_faces(ix-1,iy,iz).x();
                            sub_count(ix-1,iy,iz).x() = obs_count(ix-1,iy,iz);
                        }
                    }
                }
                break;

                case vector::Y:
                {
                    val = face_block(faceIdx).y();
                    patchId = face_patch(faceIdx).y();

                    if
                    (
                        val > limit_par
                     && iz < cellDims[vector::Z]
                     && ix < cellDims[vector::X]
                    )
                    {
                        // n_blocked_faces:
                        // count of y-faces blocked for this cell

                        if (iy < cellDims[vector::Y])
                        {
                            ++n_blocked_faces(ix,iy,iz).y();
                            sub_count(ix,iy,iz).y() = obs_count(ix,iy,iz);
                        }

                        if (iy > 0)
                        {
                            // And the neighbouring cell
                            ++n_blocked_faces(ix,iy-1,iz).y();
                            sub_count(ix,iy-1,iz).y() = obs_count(ix,iy-1,iz);
                        }
                    }
                }
                break;

                case vector::Z:
                {
                    val = face_block(faceIdx).z();
                    patchId = face_patch(faceIdx).z();

                    if
                    (
                        val > limit_par
                     && ix < cellDims[vector::X]
                     && iy < cellDims[vector::Y]
                    )
                    {
                        // n_blocked_faces:
                        // count of z-faces blocked for this cell

                        if (iz < cellDims[vector::Z])
                        {
                            ++n_blocked_faces(ix,iy,iz).z();
                            sub_count(ix,iy,iz).z() = obs_count(ix,iy,iz);
                        }

                        if (iz > 0)
                        {
                            // And the neighbouring cell
                            ++n_blocked_faces(ix,iy,iz-1).z();
                            sub_count(ix,iy,iz-1).z() = obs_count(ix,iy,iz-1);
                        }
                    }
                }
                break;
            }

            if (patchId > 0)
            {
                // If this face is on a defined patch add to list
                usedFaces(patchId).set(facei);
            }
            else if (val > limit_par)
            {
                // Add to blocked faces list
                usedFaces(PDRpatchDef::BLOCKED_FACE).set(facei);
            }
        }
    }

    // Write in time or constant dir
    const bool hasPolyMeshTimeDir = isDir(casepath/pars.timeName/"polyMesh");

    const fileName setsDir =
    (
        casepath
      / (hasPolyMeshTimeDir ? pars.timeName : word("constant"))
      / fileName("polyMesh/sets")
    );

    if (!isDir(setsDir))
    {
        mkDir(setsDir);
    }


    // Create as blockedFaces Set file for each patch, including
    // basic blocked faces
    forAll(patches, patchi)
    {
        const word& patchName = patches[patchi].patchName;

        OFstream os(setsDir / (patchName + "Set"));

        make_header(os, "polyMesh/sets", "faceSet", patchName);

        // Check for blocked faces
        const auto& fnd = usedFaces.cfind(patchi);

        if (fnd.good() && (*fnd).any())
        {
            os << nl << (*fnd).toc() << nl;
        }
        else
        {
            os << nl << labelList() << nl;
        }

        IOobject::writeEndDivider(os);
    }

    // Create the PDRMeshDict, listing the blocked faces sets and their patch names

    {
        DynamicList<word> panelNames;

        OFstream os(casepath / "system/PDRMeshDict");

        make_header(os, "system", "dictionary", "PDRMeshDict");

        os.writeEntry("blockedCells", "blockedCellsSet");
        os << nl << "blockedFaces" << nl << token::BEGIN_LIST << nl;

        for (const PDRpatchDef& p : patches)
        {
            const word& patchName = p.patchName;
            const word setName = patchName + "Set";

            if (p.patchType == 0)  // Patch
            {
                os  << "    " << token::BEGIN_LIST
                    << setName << token::SPACE
                    << patchName << token::END_LIST
                    << nl;
            }
            else if (p.patchType > 0)  // Panel
            {
                panelNames.append(setName);
            }
        }

        os  << token::END_LIST;
        os.endEntry() << nl;

        os.beginBlock("coupledFaces");

        for (const PDRpatchDef& p : patches)
        {
            const word& patchName = p.patchName;
            const word setName = patchName + "Set";

            if (p.patchType > 0)  // Panel
            {
                os.beginBlock(setName);
                os.writeEntry("wallPatch", word(patchName + "Wall"));
                os.writeEntry("cyclicMasterPatch", word(patchName + "Cyclic_half0"));
                os.endBlock();
            }
        }
        os.endBlock() << nl;

        os.writeEntry("defaultPatch", "blockedFaces");

        IOobject::writeEndDivider(os);

        // Write panelList
        OFstream os2(casepath / "panelList");

        os2<< panelNames;
        os2.endEntry();
    }
}


void write_blockedCellsSet
(
    const IjkField<scalar>& fld,
    const PDRmeshArrays& meshIndexing,
    double limit_par,
    const IjkField<Vector<direction>>& n_blocked_faces,
    const fileName& casepath,
    const word& listName
)
{
    if (listName.empty())
    {
        return;
    }

    // Write in time or constant dir
    const bool hasPolyMeshTimeDir = isDir(casepath/pars.timeName/"polyMesh");

    const fileName path =
    (
        casepath
      / (hasPolyMeshTimeDir ? pars.timeName : word("constant"))
      / fileName("polyMesh/sets")
      / listName
    );

    if (!isDir(path.path()))
    {
        mkDir(path.path());
    }

    bitSet blockedCell;

    for (label celli=0; celli < meshIndexing.nCells(); ++celli)
    {
        const labelVector& cellIdx = meshIndexing.cellIndex[celli];

        if (!isGoodIndex(cellIdx))
        {
            continue;
        }

        if (fld(cellIdx) < limit_par)
        {
            blockedCell.set(celli);
            continue;
        }

        const Vector<direction>& blocked = n_blocked_faces(cellIdx);

        const label n_bfaces = cmptSum(blocked);

        label n_bpairs = 0;

        if (n_bfaces > 1)
        {
            for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
            {
                if (blocked[cmpt] > 1) ++n_bpairs;
            }

            #if 0
            // Extra debugging
            Info<<"block " << celli << " from "
                << blocked << " -> ("
                << n_bfaces << ' ' << n_bpairs
                << ')' << nl;
            #endif
        }

        if
        (
            n_bfaces >= pars.nFacesToBlockC
         || n_bpairs >= pars.nPairsToBlockC
        )
        {
            blockedCell.set(celli);
        }
    }


    OFstream os(path);
    make_header(os, "constant/polyMesh/sets", "cellSet", listName);

    if (blockedCell.any())
    {
        os << blockedCell.toc();
    }
    else
    {
        os << labelList();
    }

    os.endEntry();

    IOobject::writeEndDivider(os);
}


// ************************************************************************* //
