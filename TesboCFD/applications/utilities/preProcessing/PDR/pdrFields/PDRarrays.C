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

#include "PDRarrays.H"
#include "PDRblock.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Smaller helper to resize matrix and assign to Zero
template<class T>
inline void resizeMatrix(SquareMatrix<T>& mat, const label n)
{
    mat.setSize(n);
    mat = Zero;
}


// Smaller helper to resize i-j-k field and assign to uniform value,
// normally Zero
template<class T>
inline void resizeField
(
    IjkField<T>& fld,
    const labelVector& ijk,
    const T& val = T(Zero)
)
{
    fld.resize(ijk);
    fld = val;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRarrays::PDRarrays()
:
    pdrBlock_(std::cref<PDRblock>(PDRblock::null()))
{}


Foam::PDRarrays::PDRarrays(const PDRblock& pdrBlock)
:
    PDRarrays()
{
    reset(pdrBlock);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRarrays::reset(const PDRblock& pdrBlock)
{
    pdrBlock_ = std::cref<PDRblock>(pdrBlock);

    // Resize all the major arrays, which are grouped in the structure arrp
    // All the relevant dimensions are in PDRblock

    // Cell-based addressing
    const labelVector cellDims = pdrBlock.sizes();

    // Face or point-based addressing
    const labelVector faceDims(cellDims + labelVector::one);

    // Max addressing dimensions for 2D arrays, with some extra space
    // These will be used for any combination of x,y,z,
    // so need to be dimensioned to the maximum size in both directions
    const label maxDim = cmptMax(pdrBlock.sizes()) + 2;

    resizeField(v_block, cellDims);
    resizeField(surf, cellDims);

    resizeField(area_block_s, cellDims);
    resizeField(area_block_r, cellDims);
    resizeField(dirn_block, cellDims);

    resizeField(face_block, faceDims);

    resizeField(along_block, cellDims);

    resizeField(betai_inv1, cellDims);

    resizeField(obs_count, cellDims);
    resizeField(sub_count, cellDims);
    resizeField(grating_count, cellDims);

    resizeField(drag_s, cellDims);
    resizeField(drag_r, cellDims);

    resizeField(obs_size, cellDims);

    for (auto& list : overlap_1d)
    {
        list.resize(maxDim);
        list = Zero;
    }

    resizeMatrix(aboverlap, maxDim);
    resizeMatrix(abperim, maxDim);
    resizeMatrix(a_lblock, maxDim);
    resizeMatrix(b_lblock, maxDim);
    resizeMatrix(ac_lblock, maxDim);
    resizeMatrix(bc_lblock, maxDim);
    resizeMatrix(c_count, maxDim);
    resizeMatrix(c_drag, maxDim);

    resizeField(face_patch, faceDims, labelVector::uniform(-1));
    resizeField(hole_in_face, faceDims);
}


// ************************************************************************* //
