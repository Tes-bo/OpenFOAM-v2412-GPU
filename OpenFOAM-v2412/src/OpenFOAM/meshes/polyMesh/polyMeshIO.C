/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::setInstance
(
    const fileName& inst,
    const IOobjectOption::writeOption wOpt
)
{
    DebugInFunction << "Resetting file instance to " << inst << endl;

    points_.writeOpt(wOpt);
    points_.instance() = inst;

    faces_.writeOpt(wOpt);
    faces_.instance() = inst;

    owner_.writeOpt(wOpt);
    owner_.instance() = inst;

    neighbour_.writeOpt(wOpt);
    neighbour_.instance() = inst;

    boundary_.writeOpt(wOpt);
    boundary_.instance() = inst;

    pointZones_.writeOpt(wOpt);
    pointZones_.instance() = inst;

    faceZones_.writeOpt(wOpt);
    faceZones_.instance() = inst;

    cellZones_.writeOpt(wOpt);
    cellZones_.instance() = inst;

    if (tetBasePtIsPtr_)
    {
        tetBasePtIsPtr_->writeOpt(wOpt);
        tetBasePtIsPtr_->instance() = inst;
    }
}


Foam::polyMesh::readUpdateState Foam::polyMesh::readUpdate()
{
    DebugInFunction << "Updating mesh based on saved data." << endl;

    // Find point/faces instances
    const fileName pointsInst(time().findInstance(meshDir(), "points"));
    const fileName facesInst(time().findInstance(meshDir(), "faces"));
    //const fileName boundInst
    //(time().findInstance(meshDir(), "boundary", IOobject::MUST_READ, facesInst));

    if (debug)
    {
        Info<< "Faces instance: old = " << facesInstance()
            << " new = " << facesInst << nl
            << "Points instance: old = " << pointsInstance()
            << " new = " << pointsInst << endl;
    }

    if (facesInst != facesInstance())
    {
        // Topological change
        if (debug)
        {
            Info<< "Topological change" << endl;
        }

        clearOut();

        // Set instance to new instance. Note that points instance can differ
        // from from faces instance.
        setInstance(facesInst);
        points_.instance() = pointsInst;

        points_.clear();
        points_ = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        faces_.clear();
        faces_ = faceCompactIOList
        (
            IOobject
            (
                "faces",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        // owner
        {
            owner_.clear();

            labelIOList list
            (
                IOobject
                (
                    "owner",
                    facesInst,
                    meshSubDir,
                    *this,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            // Update owner headerClassName.
            // The "cells" logic below may rely on it!

            owner_ = std::move(static_cast<labelList&>(list));
            owner_.headerClassName() = std::move(list.headerClassName());
            owner_.note() = std::move(list.note());
        }

        // neighbour
        {
            neighbour_.clear();

            labelIOList list
            (
                IOobject
                (
                    "neighbour",
                    facesInst,
                    meshSubDir,
                    *this,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            // Update neighbour headerClassName.
            // - not currently needed, but for symmetry with owner
            // The "cells" logic below may rely on it!

            neighbour_ = std::move(static_cast<labelList&>(list));
            neighbour_.headerClassName() = std::move(list.headerClassName());
            neighbour_.note() = std::move(list.note());
        }

        // Reset the boundary patches
        polyBoundaryMesh newBoundary
        (
            IOobject
            (
                "boundary",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this
        );

        // Check that patch types and names are unchanged
        bool boundaryChanged = false;

        if (newBoundary.size() != boundary_.size())
        {
            boundaryChanged = true;
        }
        else
        {
            forAll(boundary_, patchi)
            {
                const auto& oldPatch = boundary_[patchi];
                const auto& newPatch = newBoundary[patchi];

                if
                (
                    (oldPatch.name() != newPatch.name())
                 || (oldPatch.type() != newPatch.type())
                )
                {
                    boundaryChanged = true;
                    break;
                }
            }
        }

        if (boundaryChanged)
        {
            WarningInFunction
                << "Number of patches has changed.  This may have "
                << "unexpected consequences.  Proceed with care." << endl;

            boundary_.resize_null(newBoundary.size());

            forAll(newBoundary, patchi)
            {
                boundary_.set(patchi, newBoundary[patchi].clone(boundary_));
            }
        }
        else
        {
            forAll(boundary_, patchi)
            {
                boundary_[patchi] = polyPatch
                (
                    newBoundary[patchi].name(),
                    newBoundary[patchi].size(),
                    newBoundary[patchi].start(),
                    patchi,
                    boundary_,
                    newBoundary[patchi].physicalType(),
                    newBoundary[patchi].inGroups()
                );
            }
        }


        // Boundary is set so can use initMesh now (uses boundary_ to
        // determine internal and active faces)

        if (owner_.hasHeaderClass())
        {
            initMesh();
        }
        else
        {
            cellCompactIOList cells
            (
                IOobject
                (
                    "cells",
                    facesInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            // Recalculate the owner/neighbour addressing and reset the
            // primitiveMesh
            initMesh(cells);
        }


        // Even if number of patches stayed same still recalculate boundary
        // data.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        // Derived info
        bounds_ = boundBox(points_);
        geometricD_ = Zero;
        solutionD_ = Zero;


        // Update point/face/cell zones, but primarily just the addressing.
        // - this will be extremely fragile (not just here) if the names
        //   or the order of the zones also change

        #undef  update_meshZones
        #define update_meshZones(DataMember)                                  \
        {                                                                     \
            (DataMember).clearAddressing();                                   \
            (DataMember).clearPrimitives();                                   \
                                                                              \
            decltype(DataMember) newZones                                     \
            (                                                                 \
                IOobject                                                      \
                (                                                             \
                    (DataMember).name(),                                      \
                    facesInst,                                                \
                    meshSubDir,                                               \
                    *this,                                                    \
                    IOobject::READ_IF_PRESENT,                                \
                    IOobject::NO_WRITE,                                       \
                    IOobject::NO_REGISTER                                     \
                ),                                                            \
                *this,                                                        \
                PtrList<entry>()                                              \
            );                                                                \
            const label numZones = newZones.size();                           \
            (DataMember).resize(numZones);                                    \
                                                                              \
            for (label zonei = 0; zonei < numZones; ++zonei)                  \
            {                                                                 \
                /* Existing or new empty zone */                              \
                auto& zn = (DataMember).try_emplace                           \
                (                                                             \
                    zonei,                                                    \
                    newZones[zonei], Foam::zero{}, (DataMember)               \
                );                                                            \
                                                                              \
                /* Set addressing */                                          \
                zn.resetAddressing(std::move(newZones[zonei]));               \
            }                                                                 \
        }

        update_meshZones(pointZones_);
        update_meshZones(faceZones_);
        update_meshZones(cellZones_);
        #undef update_meshZones


        // Re-read tet base points
        tetBasePtIsPtr_ = readTetBasePtIs();


        if (boundaryChanged)
        {
            return polyMesh::TOPO_PATCH_CHANGE;
        }
        else
        {
            return polyMesh::TOPO_CHANGE;
        }
    }
    else if (pointsInst != pointsInstance())
    {
        // Points moved
        if (debug)
        {
            Info<< "Point motion" << endl;
        }

        pointIOField newPoints
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        // Re-read tet base points
        autoPtr<labelIOList> newTetBasePtIsPtr = readTetBasePtIs();

        // Update all geometry
        updateGeomPoints(std::move(newPoints), newTetBasePtIsPtr);

        return polyMesh::POINTS_MOVED;
    }
    else
    {
        if (debug)
        {
            Info<< "No change" << endl;
        }

        return polyMesh::UNCHANGED;
    }
}


// ************************************************************************* //
