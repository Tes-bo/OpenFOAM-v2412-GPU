/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

bool Foam::VF::voxel::outOfBounds
(
    const labelVector& ijk,
    const label dir
) const
{
    return (ijk[dir] < 0 || ijk[dir] >= nijk_[dir]);
};


Foam::point Foam::VF::voxel::localPosition(const vector& globalPosition) const
{
    return globalPosition - bb0_.min();
}


Foam::point Foam::VF::voxel::globalPosition(const vector& localPosition) const
{
    return bb0_.min() + localPosition;
}


void Foam::VF::voxel::setVoxelDims(const label i, const label j, const label k)
{
    nijk_[0] = max(1, i);
    nijk_[1] = max(1, j);
    nijk_[2] = max(1, k);

    dxyz_[0] = span0_[0]/nijk_[0];
    dxyz_[1] = span0_[1]/nijk_[1];
    dxyz_[2] = span0_[2]/nijk_[2];
}


void Foam::VF::voxel::refineVoxelDims()
{
    nijk_ *= 2;

    // Do not refine empty direction for 2D
    const auto& solutionD = mesh_.solutionD();
    for (direction d=0; d<3; ++d)
    {
        if (solutionD[d] == -1)
        {
            nijk_[d] = 1;
        }
    }

    setVoxelDims(nijk_[0], nijk_[1], nijk_[2]);
}


Foam::point Foam::VF::voxel::voxelMin
(
    const label i,
    const label j,
    const label k
) const
{
    return point(i*dxyz_[0], j*dxyz_[1], k*dxyz_[2]);
}


Foam::point Foam::VF::voxel::voxelMax
(
    const label i,
    const label j,
    const label k
) const
{
    return voxelMin(i+1, j+1, k+1);
}


constexpr Foam::label Foam::VF::voxel::sign0(const scalar x) const
{
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
};


Foam::labelVector Foam::VF::voxel::nijk() const noexcept
{
    return nijk_;
}


Foam::label Foam::VF::voxel::nVoxel() const noexcept
{
    return nijk_[0]*nijk_[1]*nijk_[2];
}


Foam::label Foam::VF::voxel::voxeli
(
    const labelVector ijk
) const noexcept
{
    return voxeli(ijk[0], ijk[1], ijk[2]);
}


Foam::label Foam::VF::voxel::voxeli
(
    const label i,
    const label j,
    const label k
) const noexcept
{
    return i + (nijk_[0]*(j + (nijk_[1]*k)));
}


Foam::labelVector Foam::VF::voxel::ijk(const label voxeli) const noexcept
{
    const label nx = nijk_[0];
    const label ny = nijk_[1];

    const label i = voxeli%nx;
    const label k = voxeli/nx/ny;
    const label j = (voxeli/nx)%ny;

    return labelVector(i, j, k);
}


// ************************************************************************* //
