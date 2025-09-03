/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "PDRmeshArrays.H"
#include "PDRblock.H"
#include "polyMesh.H"
#include "Time.H"
#include "IjkField.H"

// Notes
//
// Determines the face and cell numbers of all faces and cells in the
// central rectangular region where CAD_PDR operates. First,
// "points" is read and the coordinates (by which I mean here the
// indices in the x, y and z coordinate arrays) are determined. Then
// "faces" is read and for each the coordinates of the lower- x,y,z
// corner are determioned, also the orientation (X, Y or Z).
// (Orientation in the sense of e.g. + or -x is not noted.) The files
// "owner" and "neighbour" specify the six faces around each cell, so
// from these the coordinates of the cells are determined.
//
// Full checks are made that the mesh in the central region is consistent
// with CAD_PDR's mesh specified by the PDRmeshSpec file.
//
// Eventually, when writing out results, we shall work through the
// full list of cells, writing default values for any cells that are
// not in the central regtion.


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::PDRmeshArrays::gridPointRelTol = 0.02;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// A good ijk index has all components >= 0
static inline bool isGoodIndex(const Foam::labelVector& idx)
{
    return (idx.x() >= 0 && idx.y() >= 0 && idx.z() >= 0);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRmeshArrays::classify
(
    const polyMesh& mesh,
    const PDRblock& pdrBlock
)
{
    // Additional copy of i-j-k addressing
    cellDims = pdrBlock.sizes();
    faceDims = (cellDims + labelVector::one);

    const label maxPointId = cmptMax(pdrBlock.sizes())+1;

    Info<< "Mesh" << nl
        << "  nPoints:" << mesh.nPoints()
        << "  nCells:" << mesh.nCells()
        << "  nFaces:" << mesh.nFaces() << nl;

    Info<< "PDRblock" << nl
        << "  nPoints:" << pdrBlock.nPoints()
        << "  nCells:" << pdrBlock.nCells()
        << "  nFaces:" << pdrBlock.nFaces() << nl
        << "  min-edge:" << pdrBlock.edgeLimits().min() << nl;

    Info<< "Classifying ijk indexing... " << nl;


    // Bin cells into i-j-k locations with the PDRblock::findCell()
    // method, which combines a bounding box rejection and binary
    // search in the three directions.

    cellIndex.resize(mesh.nCells());
    {
        const pointField& cc = mesh.cellCentres();

        for (label celli = 0; celli < mesh.nCells(); ++celli)
        {
            cellIndex[celli] = pdrBlock.findCell(cc[celli]);
        }
    }

    // Verify that all i-j-k cells were indeed found
    {
        // This could be more efficient - but we want to be picky
        IjkField<bool> cellFound(pdrBlock.sizes(), false);

        for (label celli=0; celli < cellIndex.size(); ++celli)
        {
            const labelVector& cellIdx = cellIndex[celli];

            if (isGoodIndex(cellIdx))
            {
                cellFound(cellIdx) = true;
            }
        }

        const label firstMiss = cellFound.find(false);

        if (firstMiss >= 0)
        {
            label nMissing = 0;
            for (label celli = firstMiss; celli < cellFound.size(); ++celli)
            {
                if (!cellFound[celli])
                {
                    ++nMissing;
                }
            }

            FatalErrorInFunction
                << "No ijk location found for "
                << nMissing << " cells.\nFirst miss at: "
                << pdrBlock.index(firstMiss)
                << " indexing" << nl
                << exit(FatalError);
        }
    }


    // Bin all mesh points into i-j-k locations
    List<labelVector> pointIndex(mesh.nPoints());

    for (label pointi = 0; pointi < mesh.nPoints(); ++pointi)
    {
        const point& p = mesh.points()[pointi];
        pointIndex[pointi] = pdrBlock.gridIndex(p, gridPointRelTol);
    }

    // Min x,y,z index
    const labelMinMax invertedLimits(maxPointId, -maxPointId);
    Vector<labelMinMax> faceLimits;

    const Vector<direction> faceBits
    (
        boundBox::XDIR,
        boundBox::YDIR,
        boundBox::ZDIR
    );

    faceIndex.resize(mesh.nFaces());
    faceOrient.resize(mesh.nFaces());

    for (label facei=0; facei < mesh.nFaces(); ++facei)
    {
        labelVector& faceIdx = faceIndex[facei];

        // Check faces that are associated with i-j-k cells
        const label own = mesh.faceOwner()[facei];
        const label nei =
        (
            facei < mesh.nInternalFaces()
          ? mesh.faceNeighbour()[facei]
          : own
        );

        if (!isGoodIndex(cellIndex[own]) && !isGoodIndex(cellIndex[nei]))
        {
            // Invalid
            faceIdx.x() = faceIdx.y() = faceIdx.z() = -1;
            faceOrient[facei] = vector::X;
            continue;
        }


        faceLimits.x() = faceLimits.y() = faceLimits.z() = invertedLimits;

        for (const label pointi : mesh.faces()[facei])
        {
            for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
            {
                faceLimits[cmpt].add(pointIndex[pointi][cmpt]);
            }
        }

        direction inPlane(0u);

        for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
        {
            const auto& limits = faceLimits[cmpt];

            if (!limits.good())
            {
                // This should be impossible
                FatalErrorInFunction
                    << "Unexpected search failure for " << facei << " in "
                    << vector::componentNames[cmpt] << "-direction" << nl
                    << exit(FatalError);
            }

            if (limits.min() < 0)
            {
                FatalErrorInFunction
                    << "Face " << facei << " contains non-grid point in "
                    << vector::componentNames[cmpt] << "-direction" << nl
                    << mesh.faces()[facei] << ' '
                    << mesh.faces()[facei].points(mesh.points())
                    << exit(FatalError);
            }
            else if (limits.min() == limits.max())
            {
                // In plane
                inPlane |= faceBits[cmpt];
            }
            else if (limits.min() + 1 != limits.max())
            {
                FatalErrorInFunction
                    << "Face " << facei << " not in "
                    << vector::componentNames[cmpt] << "-plane" << nl
                    << exit(FatalError);
            }
        }

        switch (inPlane)
        {
            case boundBox::XDIR:
                faceOrient[facei] = vector::X;
                break;

            case boundBox::YDIR:
                faceOrient[facei] = vector::Y;
                break;

            case boundBox::ZDIR:
                faceOrient[facei] = vector::Z;
                break;

            default:
                FatalErrorInFunction
                    << "Face " << facei << " not in an x/y/z plane?" << nl
                    << exit(FatalError);
                break;
        }

        faceIdx.x() = faceLimits.x().min();
        faceIdx.y() = faceLimits.y().min();
        faceIdx.z() = faceLimits.z().min();
    }
}


void Foam::PDRmeshArrays::read
(
    const Time& runTime,
    const PDRblock& pdrBlock
)
{
    #include "createPolyMesh.H"
    classify(mesh, pdrBlock);
}


// ************************************************************************* //
