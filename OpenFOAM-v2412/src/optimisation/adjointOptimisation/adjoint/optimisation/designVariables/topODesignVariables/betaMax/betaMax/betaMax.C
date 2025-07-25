/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "betaMax.H"
#include "EdgeMap.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(betaMax, 0);
    defineRunTimeSelectionTable(betaMax, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::betaMax::computeLength(const dictionary& dict) const
{
    scalar length = Zero;
    // If length is not provided explicitly, loop over the provided patches
    // and compute the hydraulic diamater
    const dictionary& DarcyDict = dict.subDict(type() + "Coeffs");
    if (!DarcyDict.readIfPresent("length", length))
    {
        const labelHashSet inletPatches =
            mesh_.boundaryMesh().patchSet
            (
                DarcyDict.get<wordRes>("inletPatches")
            );

        // If 2D, use the inlet area divided by the depth in the empty direction
        if (mesh_.nGeometricD() != label(3))
        {
            // Accumulate area
            for (const label pI : inletPatches)
            {
                const fvPatch& patch = mesh_.boundary()[pI];
                length += gSum(patch.magSf());
            }

            // Divide with the span in the empty direction
            const Vector<label>& geometricD = mesh_.geometricD();
            const boundBox& bounds = mesh_.bounds();
            forAll(geometricD, iDir)
            {
                if (geometricD[iDir] == -1)
                {
                    length /= bounds.span()[iDir];
                }
            }
        }
        // If 3D, use the inlet hydraulic diameter
        else
        {
            scalar perimeter = Zero;
            scalar area = Zero;
            for (const label pI : inletPatches)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[pI];
                // Accumulate perimeter
                const edgeList& edges = patch.edges();
                const vectorField& points = patch.localPoints();
                const label nInternalEdges = patch.nInternalEdges();
                const label nEdges = patch.nEdges();
                // Processor edges should not be accounted for
                boolList isProcessorEdge = markProcessorEdges(patch);
                label nActiveEdges(0);
                forAll(isProcessorEdge, beI)
                {
                    if (!isProcessorEdge[beI])
                    {
                        perimeter += edges[nInternalEdges + beI].mag(points);
                        nActiveEdges++;
                    }
                }

                if (debug > 1)
                {
                    Info<< "patch:: " << patch.name() << nl
                        << "Number of boundary edges "
                        << returnReduce(nEdges - nInternalEdges, sumOp<label>())
                        << ", Number of non-processor edges "
                        << returnReduce(nActiveEdges, sumOp<label>()) << endl;
                }

                // Accumulate area
                area += sum(patch.magFaceAreas());
            }

            reduce(perimeter, sumOp<scalar>());
            reduce(area, sumOp<scalar>());

            length = scalar(4)*area/perimeter;
        }
    }

    return length;
}


Foam::boolList Foam::betaMax::markProcessorEdges
(
    const polyPatch& patch
) const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nNonProcessor = pbm.nNonProcessor();

    // Processor edges will artificially increase the perimeter of the inlet.
    // We need to purge them.

    // Mark all edges connected to a processor patch
    label nProcEdges(0);
    for (label procI = nNonProcessor; procI < pbm.size() ; ++procI)
    {
        const polyPatch& procPatch = pbm[procI];
        nProcEdges += procPatch.nEdges() - procPatch.nInternalEdges();
    }

    EdgeMap<bool> isInletEdge(nProcEdges);

    for (label procI = nNonProcessor; procI < pbm.size() ; ++procI)
    {
        const polyPatch& procPatch = pbm[procI];
        const labelList& procMp = procPatch.meshPoints();
        const label procInternalEdges = procPatch.nInternalEdges();
        const label procEdges = procPatch.nEdges();

        for (label edgeI = procInternalEdges; edgeI < procEdges; ++edgeI)
        {
            const edge& e = procPatch.edges()[edgeI];
            const edge meshE = edge(procMp[e[0]], procMp[e[1]]);
            isInletEdge.insert(meshE, false);
        }
    }

    // Loop over inlet edges to mark common edges with processor patches
    const label nInternalEdges = patch.nInternalEdges();
    const label nEdges = patch.nEdges();

    const edgeList& edges = patch.edges();
    const labelList& mp = patch.meshPoints();
    for (label edgeI = nInternalEdges; edgeI < nEdges; ++edgeI)
    {
        const edge& e = edges[edgeI];
        const edge meshE = edge(mp[e[0]], mp[e[1]]);
        auto iter = isInletEdge.find(meshE);

        if (iter.good())
        {
            iter.val() = true;
        }
    }

    // A special case is a processor patch intersecting the inlet patches on a
    // (true) boundary edge of the latter. Use syncEdgeMap to make sure these
    // edges don't make it to the final list
    syncTools::syncEdgeMap
    (
        pbm.mesh(),
        isInletEdge,
        andEqOp<bool>()
    );

    // Now loop over the inlet faces once more and mark the common edges with
    // processor boundaries
    boolList isProcessorEdge(nEdges - nInternalEdges, false);
    for (label edgeI = nInternalEdges; edgeI < nEdges; ++edgeI)
    {
        const edge& e = edges[edgeI];
        const edge meshE = edge(mp[e[0]], mp[e[1]]);

        if (isInletEdge.lookup(meshE, false))
        {
            isProcessorEdge[edgeI - nInternalEdges] = true;
        }
    }

    return isProcessorEdge;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::betaMax::betaMax
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    value_(dict.getOrDefault<scalar>("betaMax", Zero))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::betaMax> Foam::betaMax::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.getOrDefault<word>("betaMaxType", "value"));

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    Info<< "betaMax type " << modelType << endl;

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "betaMaxType",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<betaMax>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::betaMax::value() const
{
    return value_;
}


// ************************************************************************* //
