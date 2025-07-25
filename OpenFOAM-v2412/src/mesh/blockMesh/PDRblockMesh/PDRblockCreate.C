/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "PDRblock.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRblock::createPoints(pointField& pts) const
{
    const label ni = sizes().x();
    const label nj = sizes().y();
    const label nk = sizes().z();

    pts.resize(nPoints());

    for (label k=0; k<=nk; ++k)
    {
        for (label j=0; j<=nj; ++j)
        {
            for (label i=0; i<=ni; ++i)
            {
                point& pt = pts[pointLabel(i,j,k)];

                pt.x() = grid_.x()[i];
                pt.y() = grid_.y()[j];
                pt.z() = grid_.z()[k];
            }
        }
    }
}


Foam::label Foam::PDRblock::addInternalFaces
(
    faceList::iterator& faceIter,
    labelList::iterator& ownIter,
    labelList::iterator& neiIter
) const
{
    const label ni = sizes().x();
    const label nj = sizes().y();
    const label nk = sizes().z();

    const labelList::iterator firstIter = ownIter;

    for (label k=0; k<nk; ++k)
    {
        for (label j=0; j<nj; ++j)
        {
            for (label i=0; i<ni; ++i)
            {
                const label celli = cellLabel(i, j, k);

                // Local Face 1 == x-max
                if (i < ni-1)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i+1, j,   k);
                    f[1] = pointLabel(i+1, j+1, k);
                    f[2] = pointLabel(i+1, j+1, k+1);
                    f[3] = pointLabel(i+1, j,   k+1);

                    *ownIter = celli;
                    *neiIter = cellLabel(i+1, j, k);

                    ++ownIter;
                    ++neiIter;
                }

                // Local Face 3 == y-max
                if (j < nj-1)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   j+1, k);
                    f[1] = pointLabel(i,   j+1, k+1);
                    f[2] = pointLabel(i+1, j+1, k+1);
                    f[3] = pointLabel(i+1, j+1, k);

                    *ownIter = celli;
                    *neiIter = cellLabel(i, j+1, k);

                    ++ownIter;
                    ++neiIter;
                }

                // Local Face 5 == z-max
                if (k < nk-1)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   j,   k+1);
                    f[1] = pointLabel(i+1, j,   k+1);
                    f[2] = pointLabel(i+1, j+1, k+1);
                    f[3] = pointLabel(i,   j+1, k+1);

                    *ownIter = celli;
                    *neiIter = cellLabel(i, j, k+1);

                    ++ownIter;
                    ++neiIter;
                }
            }
        }
    }

    // Return the number of faces added
    return (ownIter - firstIter);
}


Foam::label Foam::PDRblock::addBoundaryFaces
(
    const direction shapeFacei,
    faceList::iterator& faceIter,
    labelList::iterator& ownIter
) const
{
    const label ni = sizes().x();
    const label nj = sizes().y();
    const label nk = sizes().z();

    const labelList::iterator firstIter = ownIter;

    switch (shapeFacei)
    {
        // Face 0 == x-min
        case 0:
        {
            for (label k=0; k<nk; ++k)
            {
                for (label j=0; j<nj; ++j)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(0, j,   k);
                    f[1] = pointLabel(0, j,   k+1);
                    f[2] = pointLabel(0, j+1, k+1);
                    f[3] = pointLabel(0, j+1, k);

                    *ownIter = cellLabel(0, j, k);
                    ++ownIter;
                }
            }
        }
        break;

        // Face 1 == x-max
        case 1:
        {
            for (label k=0; k<nk; ++k)
            {
                for (label j=0; j<nj; ++j)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(ni, j,   k);
                    f[1] = pointLabel(ni, j+1, k);
                    f[2] = pointLabel(ni, j+1, k+1);
                    f[3] = pointLabel(ni, j,   k+1);

                    *ownIter = cellLabel(ni-1, j, k);
                    ++ownIter;
                }
            }
        }
        break;

        // Face 2 == y-min
        case 2:
        {
            for (label i=0; i<ni; ++i)
            {
                for (label k=0; k<nk; ++k)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   0, k);
                    f[1] = pointLabel(i+1, 0, k);
                    f[2] = pointLabel(i+1, 0, k+1);
                    f[3] = pointLabel(i,   0, k+1);

                    *ownIter = cellLabel(i, 0, k);
                    ++ownIter;
                }
            }
        }
        break;

        // Face 3 == y-max
        case 3:
        {
            for (label i=0; i<ni; ++i)
            {
                for (label k=0; k<nk; ++k)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   nj, k);
                    f[1] = pointLabel(i,   nj, k+1);
                    f[2] = pointLabel(i+1, nj, k+1);
                    f[3] = pointLabel(i+1, nj, k);

                    *ownIter = cellLabel(i, nj-1, k);
                    ++ownIter;
                }
            }
        }
        break;

        // Face 4 == z-min
        case 4:
        {
            for (label i=0; i<ni; ++i)
            {
                for (label j=0; j<nj; ++j)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   j,   0);
                    f[1] = pointLabel(i,   j+1, 0);
                    f[2] = pointLabel(i+1, j+1, 0);
                    f[3] = pointLabel(i+1, j,   0);

                    *ownIter = cellLabel(i, j, 0);
                    ++ownIter;
                }
            }
        }
        break;

        // Face 5 == z-max
        case 5:
        {
            for (label i=0; i<ni; ++i)
            {
                for (label j=0; j<nj; ++j)
                {
                    auto& f = *faceIter;
                    ++faceIter;
                    f.resize(4);

                    f[0] = pointLabel(i,   j,   nk);
                    f[1] = pointLabel(i+1, j,   nk);
                    f[2] = pointLabel(i+1, j+1, nk);
                    f[3] = pointLabel(i,   j+1, nk);

                    *ownIter = cellLabel(i, j, nk-1);
                    ++ownIter;
                }
            }
        }
        break;
    }

    // Return the number of faces added
    return (ownIter - firstIter);
}


Foam::autoPtr<Foam::polyMesh>
Foam::PDRblock::innerMesh(const IOobject& io) const
{
    pointField pts(nPoints());

    faceList faces(nFaces());
    labelList own(nFaces());
    labelList nei(nInternalFaces());

    auto faceIter = faces.begin();
    auto ownIter  = own.begin();
    auto neiIter  = nei.begin();

    createPoints(pts);

    addInternalFaces(faceIter, ownIter, neiIter);

    // After readBoundary() we have a complete list of patches
    // without any conflicts, and the correct size per-patch

    // Add boundary faces and adjust patch sizes
    for (const boundaryEntry& bentry : patches_)
    {
        for (const label shapeFacei : bentry.faces_)
        {
            addBoundaryFaces(shapeFacei, faceIter, ownIter);
        }
    }


    IOobject iomesh(io);
    iomesh.writeOpt(IOobject::AUTO_WRITE);

    auto meshPtr = autoPtr<polyMesh>::New
    (
        iomesh,
        std::move(pts),
        std::move(faces),
        std::move(own),
        std::move(nei)
    );
    polyMesh& pmesh = *meshPtr;

    polyPatchList patches(patches_.size());

    label startFace = nInternalFaces();

    label patchi = 0;

    for (const boundaryEntry& bentry : patches_)
    {
        patches.set
        (
            patchi,
            polyPatch::New
            (
                bentry.type_,
                bentry.name_,
                bentry.size_,
                startFace,
                patchi,  // index
                pmesh.boundaryMesh()
            )
        );

        // physicalType?

        startFace += bentry.size_;
        ++patchi;
    }

    pmesh.addPatches(patches);

    return meshPtr;
}


Foam::autoPtr<Foam::polyMesh>
Foam::PDRblock::mesh(const IOobject& io) const
{
    if (outer_.active())
    {
        Info<< "Outer region is active, using blockMesh generation" << nl;
        return meshBlockMesh(io);
    }
    else
    {
        Info<< "Outer region is inactive, using ijk generation" << nl;
        return innerMesh(io);
    }
}


// ************************************************************************* //
