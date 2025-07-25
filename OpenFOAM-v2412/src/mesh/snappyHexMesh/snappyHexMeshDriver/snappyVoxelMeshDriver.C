/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "snappyVoxelMeshDriver.H"
#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "refinementParameters.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "searchableSurfaces.H"
#include "voxelMeshSearch.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(snappyVoxelMeshDriver, 0);
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::snappyVoxelMeshDriver::addNeighbours
(
    const labelList& cellLevel,
    const labelVector& voxel,
    const label voxeli,
    DynamicList<labelVector>& front
) const
{
    const labelVector off(voxelMeshSearch::offset(n_));

    for (direction dir = 0; dir < 3; ++dir)
    {
        if (voxel[dir] > 0)
        {
            labelVector left(voxel);
            left[dir] -= 1;
            if (cellLevel[voxeli-off[dir]] == -1)
            {
                front.append(left);
            }
        }
        if (voxel[dir] < n_[dir]-1)
        {
            labelVector right(voxel);
            right[dir] += 1;
            if (cellLevel[voxeli+off[dir]] == -1)
            {
                front.append(right);
            }
        }
    }
}


// Insert cell level for all volume refinement
Foam::tmp<Foam::pointField> Foam::snappyVoxelMeshDriver::voxelCentres() const
{
    tmp<pointField> tcc(tmp<pointField>::New(n_.x()*n_.y()*n_.z()));
    pointField& cc = tcc.ref();

    const labelVector off(voxelMeshSearch::offset(n_));
    label voxeli = voxelMeshSearch::index(n_, labelVector(0, 0, 0));
    for (label k = 0; k < n_[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = 0; j < n_[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = 0; i < n_[0]; i++)
            {
                const labelVector voxel(i, j, k);
                cc[voxeli] = voxelMeshSearch::centre(bb_, n_, voxel);
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
    return tcc;
}


void Foam::snappyVoxelMeshDriver::isInside
(
    const pointField& cc,
    boolList& isVoxelInMesh
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    isVoxelInMesh.setSize(cc.size());
    if (isVoxelInMesh.size() < mesh.globalData().nTotalCells())
    {
        forAll(cc, voxeli)
        {
            const label celli = mesh.findCell
            (
                cc[voxeli],
                polyMesh::FACE_PLANES
            );
            isVoxelInMesh[voxeli] = (celli != -1);
        }
        Pstream::listCombineGather(isVoxelInMesh, orEqOp<bool>());
    }
    else
    {
        for (label celli = 0; celli < mesh.nCells(); ++celli)
        {
            const boundBox cellBb(mesh.cellBb(celli));

            voxelMeshSearch::fill
            (
                isVoxelInMesh,
                bb_,
                n_,
                cellBb,
                1,
                orEqOp<bool>()
            );
        }
        Pstream::listCombineGather(isVoxelInMesh, orEqOp<bool>());
    }
}


void Foam::snappyVoxelMeshDriver::markSurfaceRefinement
(
    labelList& voxelLevel,
    labelList& globalRegion
) const
{
    // Insert cell level for all refinementSurfaces

    const refinementSurfaces& s = meshRefiner_.surfaces();
    forAll(s.surfaces(), surfi)
    {
        label geomi = s.surfaces()[surfi];
        const searchableSurface& geom = s.geometry()[geomi];
        //Pout<< "Geometry:" << s.names()[surfi] << endl;
        if (isA<triSurface>(geom))
        {
            const triSurface& ts = refCast<const triSurface>(geom);
            const pointField& points = ts.points();

            for (const labelledTri& tri : ts)
            {
                label regioni = tri.region();
                label globalRegioni = s.regionOffset()[surfi]+regioni;
                const boundBox triBb(tri.box(points));

                // Fill cellLevel
                label level = s.minLevel()[globalRegioni];
                voxelMeshSearch::fill
                (
                    voxelLevel,
                    bb_,
                    n_,
                    triBb,
                    level,
                    maxEqOp<label>()
                );
                voxelMeshSearch::fill
                (
                    globalRegion,
                    bb_,
                    n_,
                    triBb,
                    globalRegioni,
                    maxEqOp<label>()
                );
            }
        }
        // else: maybe do intersection tests?
    }
}


void Foam::snappyVoxelMeshDriver::findVoxels
(
    const labelList& voxelLevel,
    const pointField& locationsOutsideMesh,
    labelList& voxels
) const
{
    voxels.setSize(locationsOutsideMesh.size());
    voxels = -1;
    forAll(locationsOutsideMesh, loci)
    {
        const point& pt = locationsOutsideMesh[loci];
        label voxeli = voxelMeshSearch::index(bb_, n_, pt, false);

        if (voxeli == -1 || voxelLevel[voxeli] == labelMax)
        {
            WarningInFunction << "Location outside mesh "
                << pt << " is outside mesh with bounding box "
                << bb_ << endl;
        }
        else
        {
            voxels[loci] = voxeli;
        }
    }
}


void Foam::snappyVoxelMeshDriver::floodFill
(
    const label startVoxeli,
    const label newLevel,
    labelList& voxelLevel
) const
{
    DynamicList<labelVector> front;
    front.append(voxelMeshSearch::index3(n_, startVoxeli));

    DynamicList<labelVector> newFront;
    while (true)
    {
        newFront.clear();
        for (const auto& voxel : front)
        {
            label voxeli = voxelMeshSearch::index(n_, voxel);
            if (voxelLevel[voxeli] == -1)
            {
                voxelLevel[voxeli] = 0;
                addNeighbours
                (
                    voxelLevel,
                    voxel,
                    voxeli,
                    newFront
                );
            }
        }

        if (newFront.empty())
        {
            break;
        }
        front.transfer(newFront);
    }
}


void Foam::snappyVoxelMeshDriver::max
(
    const labelList& maxLevel,
    labelList& voxelLevel
) const
{
    // Mark voxels with level
    const labelVector off(voxelMeshSearch::offset(n_));

    label voxeli = voxelMeshSearch::index(n_, labelVector(0, 0, 0));
    for (label k = 0; k < n_[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = 0; j < n_[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = 0; i < n_[0]; i++)
            {
                voxelLevel[voxeli] = Foam::max
                (
                    voxelLevel[voxeli],
                    maxLevel[voxeli]
                );
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }
}


Foam::labelList Foam::snappyVoxelMeshDriver::count
(
    const labelList& voxelLevel
) const
{

    label maxLevel = 0;
    for (const auto level : voxelLevel)
    {
        if (level != labelMax)
        {
            maxLevel = Foam::max(maxLevel, level);
        }
    }
    labelList count(maxLevel+1, 0);

    const labelVector off(voxelMeshSearch::offset(n_));

    label voxeli = voxelMeshSearch::index(n_, labelVector(0, 0, 0));
    for (label k = 0; k < n_[2]; k++)
    {
        const label start1 = voxeli;
        for (label j = 0; j < n_[1]; j++)
        {
            const label start0 = voxeli;
            for (label i = 0; i < n_[0]; i++)
            {
                label level = voxelLevel[voxeli];

                if (level != -1 && level != labelMax)
                {
                    ++count[level];
                }
                voxeli += off[0];
            }
            voxeli = start0 + off[1];
        }
        voxeli = start1 + off[2];
    }

    return count;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyVoxelMeshDriver::snappyVoxelMeshDriver
(
    meshRefinement& meshRefiner,
    const labelUList& globalToMasterPatch,
    const labelUList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch),
    bb_(meshRefiner_.mesh().bounds())
{
    label maxLevel = labelMin;

    // Feature refinement
    const labelListList& featLevels = meshRefiner_.features().levels();
    forAll(featLevels, feati)
    {
        maxLevel = Foam::max(maxLevel, Foam::max(featLevels[feati]));
    }

    // Surface refinement
    const labelList& surfaceLevels = meshRefiner_.surfaces().maxLevel();
    maxLevel = Foam::max(maxLevel, Foam::max(surfaceLevels));

    // Shell refinement
    maxLevel = Foam::max(maxLevel, meshRefiner_.shells().maxLevel());

    const scalar level0Len = meshRefiner_.meshCutter().level0EdgeLength();

    const int oldWidth = Sout.width();

    Info<< nl
        << "Cell size estimate :" << nl
        << "    Level "
        << setw(2) << label(0) << setw(oldWidth)
        << " : " << level0Len << nl
        << "    Level "
        << setw(2) << maxLevel << setw(oldWidth)
        << " : " << level0Len/pow(2.0, maxLevel) << nl
        << endl;


    // Define voxel mesh with similar dimensions as mesh
    const vector meshSpan(bb_.span());
    n_ = labelVector
    (
        round(meshSpan.x()/level0Len),
        round(meshSpan.y()/level0Len),
        round(meshSpan.z()/level0Len)
    );
    label nTot = n_.x()*n_.y()*n_.z();
    while (nTot < 1000000)  //1048576)
    {
        n_ *= 2;
        nTot = n_.x()*n_.y()*n_.z();
    }

    Info<< "Voxellating initial mesh : " << n_ << nl << endl;

    tmp<pointField> tcc(voxelCentres());
    const pointField& cc = tcc();

    Info<< "Voxel refinement :" << nl
        << "    Initial                      : (" << nTot << ')' << endl;

    boolList isVoxelInMesh;
    isInside(cc, isVoxelInMesh);

    if (Pstream::master())
    {
        voxelLevel_.setSize(nTot, -1);
        globalRegion_.setSize(nTot, -1);

        // Remove cells outside initial mesh
        forAll(isVoxelInMesh, voxeli)
        {
            if (!isVoxelInMesh[voxeli])
            {
                voxelLevel_[voxeli] = labelMax;
                globalRegion_[voxeli] = -1;
            }
        }

        //if (debug)
        {
            Info<< "    After removing outside cells : " << count(voxelLevel_)
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::snappyVoxelMeshDriver::doRefine
(
    const refinementParameters& refineParams
)
{
    const scalar level0Len = meshRefiner_.meshCutter().level0EdgeLength();

    tmp<pointField> tcc(voxelCentres());
    const pointField& cc = tcc();

    boolList isVoxelInMesh;
    isInside(cc, isVoxelInMesh);

    if (Pstream::master())
    {
        // Mark voxels containing (parts of) triangles
        markSurfaceRefinement(voxelLevel_, globalRegion_);

        //if (debug)
        {
            Info<< "    After surface refinement     : " << count(voxelLevel_)
                << endl;
        }


        // Find outside locations (and their current refinement level)
        const pointField& outsidePoints = refineParams.locationsOutsideMesh();
        labelList outsideMeshVoxels;
        findVoxels
        (
            voxelLevel_,
            outsidePoints,
            outsideMeshVoxels
        );
        labelList outsideOldLevel(outsideMeshVoxels.size(), -1);
        forAll(outsideMeshVoxels, loci)
        {
            label voxeli = outsideMeshVoxels[loci];
            if (voxeli >= 0)
            {
                outsideOldLevel[loci] = voxelLevel_[outsideMeshVoxels[loci]];
                if (outsideOldLevel[loci] >= 0)
                {
                    WarningInFunction << "Location outside mesh "
                        << outsidePoints[loci]
                        << " is inside mesh or close to surface" << endl;
                }
            }
        }


        // Find inside locations
        labelList insideMeshVoxels;
        findVoxels
        (
            voxelLevel_,
            refineParams.locationsInMesh(),
            insideMeshVoxels
        );

        forAll(insideMeshVoxels, loci)
        {
            label voxeli = insideMeshVoxels[loci];
            if (voxeli != -1)
            {
                if (voxelLevel_[voxeli] != -1)
                {
                    WarningInFunction << "Location inside mesh "
                        << refineParams.locationsInMesh()[loci]
                        << " is marked as a surface voxel " << voxeli
                        << " with cell level " << voxelLevel_[voxeli] << endl;
                }
                else
                {
                    // Flood-fill out from voxel
                    floodFill(voxeli, 0, voxelLevel_);
                }
            }
        }

        //if (debug)
        {
            Info<< "    After keeping inside voxels  : " << count(voxelLevel_)
                << endl;
        }


        // Re-check the outside locations to see if they have been bled into
        {
            forAll(outsideMeshVoxels, loci)
            {
                label voxeli = outsideMeshVoxels[loci];
                if (voxeli >= 0 && voxelLevel_[voxeli] != outsideOldLevel[loci])
                {
                    WarningInFunction << "Location outside mesh "
                        << outsidePoints[loci]
                        << " is reachable from an inside location" << nl
                        << "Either your locations are too close to the"
                        << " geometry or there might be a leak in the"
                        << " geometry" << endl;
                }
            }
        }


        // Shell refinement : find ccs inside higher shells
        labelList maxLevel;
        meshRefiner_.shells().findHigherLevel(cc, voxelLevel_, maxLevel);

        // Assign max of maxLevel and voxelLevel
        max(maxLevel, voxelLevel_);

        // Determine number of levels
        const labelList levelCounts(count(voxelLevel_));

        //if (debug)
        {
            Info<< "    After shell refinement       : " << levelCounts << endl;
        }


        const vector meshSpan(bb_.span());
        const vector voxel0Size
        (
            meshSpan[0]/n_[0],
            meshSpan[1]/n_[1],
            meshSpan[2]/n_[2]
        );
        label cellCount = 0;
        forAll(levelCounts, leveli)
        {
            const scalar s = level0Len/pow(2.0, leveli);
            const scalar nCellsPerVoxel
            (
                voxel0Size[0]/s
               *voxel0Size[1]/s
               *voxel0Size[2]/s
            );
            cellCount += levelCounts[leveli]*nCellsPerVoxel;
        }
        Info<< "Estimated cell count : " << cellCount << endl;
    }
}


// ************************************************************************* //
