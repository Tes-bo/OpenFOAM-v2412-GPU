/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Handle merging of patch pairs

\*---------------------------------------------------------------------------*/

{
    wordPairList mergePatchPairs;

    // Read in a list of merge patch pairs
    if
    (
        meshDict.readIfPresent("mergePatchPairs", mergePatchPairs)
     && mergePatchPairs.size()
    )
    {
        Info<< "Merging " << mergePatchPairs.size() << " patch pairs" << nl;

        // Cleanup
        wordHashSet cleanupPatches(4*mergePatchPairs.size());
        wordHashSet cleanupPointZones(2*mergePatchPairs.size());
        wordHashSet cleanupFaceZones(2*mergePatchPairs.size());

        Info<< "    Adding point and face zones" << endl;
        {
            const auto& pbm = mesh.boundaryMesh();

            auto& pzs = mesh.pointZones(); pzs.clearAddressing();
            auto& fzs = mesh.faceZones();  fzs.clearAddressing();

            forAll(mergePatchPairs, pairi)
            {
                // Patch pairs
                const polyPatch& patch0 = pbm[mergePatchPairs[pairi].first()];
                const polyPatch& patch1 = pbm[mergePatchPairs[pairi].second()];

                const word mergeName
                (
                    mergePatchPairs[pairi].first()
                  + mergePatchPairs[pairi].second()
                  + Foam::name(pairi)
                );

                // An empty zone for cut points
                pzs.emplace_back
                (
                    mergeName + "CutPointZone",
                    pzs.size(),     // index
                    pzs
                );
                cleanupPointZones.insert(pzs.back().name());

                // Coupling side 0 (master)
                fzs.emplace_back
                (
                    mergeName + "Side0Zone",
                    identity(patch0.range()),
                    false,          // none are flipped
                    fzs.size(),     // index
                    fzs
                );
                cleanupFaceZones.insert(fzs.back().name());

                // Coupling side 1 (slave)
                fzs.emplace_back
                (
                    mergeName + "Side1Zone",
                    identity(patch1.range()),
                    false,          // none are flipped
                    fzs.size(),     // index
                    fzs
                );
                cleanupFaceZones.insert(fzs.back().name());

                // An empty zone for cut faces
                fzs.emplace_back
                (
                    mergeName + "CutFaceZone",
                    fzs.size(),     // index
                    fzs
                );
                cleanupFaceZones.insert(fzs.back().name());
            }
        }


        Info<< "    Merging with attachPolyTopoChanger" << endl;
        attachPolyTopoChanger polyMeshAttacher(mesh);
        polyMeshAttacher.resize(1);

        forAll(mergePatchPairs, pairi)
        {
            cleanupPatches.insert(mergePatchPairs[pairi].first());
            cleanupPatches.insert(mergePatchPairs[pairi].second());

            const word mergeName
            (
                mergePatchPairs[pairi].first()
              + mergePatchPairs[pairi].second()
              + Foam::name(pairi)
            );

            // Add the sliding interface mesh modifier
            polyMeshAttacher.set
            (
                0,
                new slidingInterface
                (
                    "couple" + Foam::name(pairi),
                    pairi,
                    polyMeshAttacher,
                    mergeName + "Side0Zone",
                    mergeName + "Side1Zone",
                    mergeName + "CutPointZone",
                    mergeName + "CutFaceZone",
                    mergePatchPairs[pairi].first(),
                    mergePatchPairs[pairi].second(),
                    slidingInterface::INTEGRAL, // always integral
                    false,
                    intersection::VISIBLE
                )
            );

            polyMeshAttacher.attach(false);  // Do not yet remove empty patches
        }

        // Re-do the boundary patches, removing empty merge patches
        // but keeping any other empty patches
        {
            const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

            polyPatchList newPatches(oldPatches.size());
            label nNewPatches = 0;

            wordHashSet removedPatches(cleanupPatches.capacity());

            forAll(oldPatches, patchi)
            {
                const word& patchName = oldPatches[patchi].name();

                if
                (
                    !cleanupPatches.found(patchName)
                 || returnReduceOr(oldPatches[patchi].size())
                )
                {
                    newPatches.set
                    (
                        nNewPatches,
                        oldPatches[patchi].clone
                        (
                            mesh.boundaryMesh(),
                            nNewPatches,
                            oldPatches[patchi].size(),
                            oldPatches[patchi].start()
                        )
                    );

                    ++nNewPatches;
                }
                else
                {
                    removedPatches.insert(patchName);
                }
            }

            newPatches.resize(nNewPatches);

            mesh.removeBoundary();
            mesh.addPatches(newPatches);

            Info<< "Removed " << removedPatches.size()
                << " empty merged patches:" << nl
                << "    " << flatOutput(removedPatches.sortedToc()) << endl;
        }

        // Cleanup empty merged point zones
        {
            PtrList<pointZone>& zones = mesh.pointZones();
            mesh.pointZones().clearAddressing();

            wordHashSet removedZones(2*zones.size());

            label nZones = 0;
            forAll(zones, zonei)
            {
                if
                (
                    !cleanupPointZones.found(zones[zonei].name())
                 || returnReduceOr(zones[zonei].size())
                )
                {
                    zones.set(nZones, zones.release(zonei));
                    zones[nZones].index() = nZones;  // re-index
                    ++nZones;
                }
                else
                {
                    removedZones.insert(zones[zonei].name());
                }
            }
            zones.resize(nZones);

            if (removedZones.size())
            {
                Info<< "Removed " << removedZones.size()
                    << " empty point zones:" << nl
                    << "    " << flatOutput(removedZones.sortedToc()) << endl;
            }
        }

        // Cleanup empty merged face zones
        {
            PtrList<faceZone>& zones = mesh.faceZones();
            mesh.faceZones().clearAddressing();

            wordHashSet removedZones(2*zones.size());

            label nZones = 0;
            forAll(zones, zonei)
            {
                if
                (
                    !cleanupFaceZones.found(zones[zonei].name())
                 || returnReduceOr(zones[zonei].size())
                )
                {
                    zones.set(nZones, zones.release(zonei));
                    zones[nZones].index() = nZones;  // re-index
                    ++nZones;
                }
                else
                {
                    removedZones.insert(zones[zonei].name());
                }
            }
            zones.resize(nZones);

            if (removedZones.size())
            {
                Info<< "Removed " << removedZones.size()
                    << " empty merged face zones:" << nl
                    << "    " << flatOutput(removedZones.sortedToc()) << endl;
            }
        }
    }
    else
    {
        Info<< "No patch pairs to merge" << endl;
    }
}


// ************************************************************************* //
