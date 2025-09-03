/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "MapLagrangianFields.H"
#include "passiveParticleCloud.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

static const scalar perturbFactor = 1e-6;


// Special version of findCell that generates a cell guaranteed to be
// compatible with tracking.
static label findCell(const Cloud<passiveParticle>& cloud, const point& pt)
{
    label celli = -1;
    label tetFacei = -1;
    label tetPtI = -1;

    const polyMesh& mesh = cloud.pMesh();

    mesh.findCellFacePt(pt, celli, tetFacei, tetPtI);

    if (celli >= 0)
    {
        return celli;
    }
    else
    {
        // See if particle on face by finding nearest face and shifting
        // particle.

        meshSearch meshSearcher
        (
            mesh,
            polyMesh::FACE_PLANES    // no decomposition needed
        );

        label facei = meshSearcher.findNearestBoundaryFace(pt);

        if (facei >= 0)
        {
            const point& cc = mesh.cellCentres()[mesh.faceOwner()[facei]];

            const point perturbPt = (1-perturbFactor)*pt+perturbFactor*cc;

            mesh.findCellFacePt(perturbPt, celli, tetFacei, tetPtI);

            return celli;
        }
    }

    return -1;
}


void mapLagrangian(const meshToMesh& interp)
{
    // Determine which particles are in meshTarget
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyMesh& meshSource = interp.srcRegion();
    const polyMesh& meshTarget = interp.tgtRegion();
    const labelListList& sourceToTarget = interp.srcToTgtCellAddr();

    const fileNameList cloudDirs
    (
        readDir
        (
            meshSource.time().timePath()/cloud::prefix,
            fileName::DIRECTORY
        )
    );

    for (const fileName& cloudDir : cloudDirs)
    {
        // Search for list of lagrangian objects for this time
        IOobjectList objects
        (
            meshSource,
            meshSource.time().timeName(),
            cloud::prefix/cloudDir
        );

        if
        (
            returnReduceOr
            (
                objects.found("coordinates") || objects.found("positions")
            )
        )
        {
            // Has coordinates/positions - so must be a valid cloud
            Info<< nl << "    processing cloud " << cloudDir << endl;

            // Read positions & cell
            passiveParticleCloud sourceParcels
            (
                meshSource,
                cloudDir,
                false
            );

            Info<< "    read " << sourceParcels.size()
                << " parcels from source mesh." << endl;

            // Construct empty target cloud
            passiveParticleCloud targetParcels
            (
                meshTarget,
                Foam::zero{},
                cloudDir
            );

            passiveParticle::trackingData td(targetParcels);

            label sourceParticleI = 0;

            // Indices of source particles that get added to targetParcels
            DynamicList<label> addParticles(sourceParcels.size());

            // Unmapped particles
            labelHashSet unmappedSource(sourceParcels.size());


            // Initial: track from fine-mesh cell centre to particle position
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // This requires there to be no boundary in the way.


            for (const passiveParticle& p : sourceParcels)
            {
                bool foundCell = false;

                // Assume that cell from read parcel is the correct one...
                if (p.cell() >= 0)
                {
                    const labelList& targetCells =
                        sourceToTarget[p.cell()];

                    // Particle probably in one of the targetcells. Try
                    // all by tracking from their cell centre to the parcel
                    // position.

                    for (const label targetCell : targetCells)
                    {
                        // Track from its cellcentre to position to make sure.
                        autoPtr<passiveParticle> newPtr
                        (
                            new passiveParticle
                            (
                                meshTarget,
                                barycentric(1, 0, 0, 0),
                                targetCell,
                                meshTarget.cells()[targetCell][0],
                                1
                            )
                        );
                        passiveParticle& newP = newPtr();

                        newP.track(p.position() - newP.position(), 0);

                        if (!newP.onFace())
                        {
                            // Hit position.
                            foundCell = true;
                            addParticles.append(sourceParticleI);
                            targetParcels.addParticle(newPtr.ptr());
                            break;
                        }
                    }
                }

                if (!foundCell)
                {
                    // Store for closer analysis
                    unmappedSource.insert(sourceParticleI);
                }

                sourceParticleI++;
            }

            Info<< "    after meshToMesh addressing found "
                << targetParcels.size()
                << " parcels in target mesh." << endl;


            // Do closer inspection for unmapped particles
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (unmappedSource.size())
            {
                sourceParticleI = 0;

                for (passiveParticle& p : sourceParcels)
                {
                    if (unmappedSource.found(sourceParticleI))
                    {
                        const label targetCell =
                            findCell(targetParcels, p.position());

                        if (targetCell >= 0)
                        {
                            unmappedSource.erase(sourceParticleI);
                            addParticles.append(sourceParticleI);
                            targetParcels.addParticle
                            (
                                new passiveParticle
                                (
                                    meshTarget,
                                    p.position(),
                                    targetCell
                                )
                            );
                            sourceParcels.remove(&p);
                        }
                    }
                    sourceParticleI++;
                }
            }
            addParticles.shrink();

            Info<< "    after additional mesh searching found "
                << targetParcels.size() << " parcels in target mesh." << endl;

            if (addParticles.size())
            {
                IOPosition<passiveParticleCloud>(targetParcels).write();

                // addParticles now contains the indices of the sourceMesh
                // particles that were appended to the target mesh.

                // Map lagrangian fields
                // ~~~~~~~~~~~~~~~~~~~~~

                MapLagrangianFields<label>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
                MapLagrangianFields<scalar>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
                MapLagrangianFields<vector>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
                MapLagrangianFields<sphericalTensor>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
                MapLagrangianFields<symmTensor>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
                MapLagrangianFields<tensor>
                (
                    cloudDir,
                    objects,
                    meshTarget,
                    addParticles
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
