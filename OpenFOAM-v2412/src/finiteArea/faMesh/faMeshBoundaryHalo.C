/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "faMeshBoundaryHalo.H"
#include "faMesh.H"
#include "globalIndex.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMeshBoundaryHalo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshBoundaryHalo::faMeshBoundaryHalo(const label comm)
:
    mapDistributeBase(comm),
    inputMeshFaces_(),
    boundaryToCompact_()
{}


Foam::faMeshBoundaryHalo::faMeshBoundaryHalo(const faMesh& areaMesh)
:
    mapDistributeBase(),
    inputMeshFaces_(),
    boundaryToCompact_()
{
    reset(areaMesh);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshBoundaryHalo::clear()
{
    static_cast<mapDistributeBase&>(*this) = mapDistributeBase();

    inputMeshFaces_.clear();
    boundaryToCompact_.clear();
}


Foam::label Foam::faMeshBoundaryHalo::haloSize() const
{
    if (Pstream::parRun())
    {
        return boundaryToCompact_.size();
    }
    else
    {
        return inputMeshFaces_.size();
    }
}


void Foam::faMeshBoundaryHalo::reset(const faMesh& areaMesh)
{
    inputMeshFaces_.clear();
    boundaryToCompact_.clear();

    const auto& procConnections = areaMesh.boundaryConnections();

    if (!Pstream::parRun())
    {
        // Serial - extract halo numbers directly

        inputMeshFaces_.resize(procConnections.size());

        forAll(procConnections, connecti)
        {
            // Connected neighbour, non-parallel = must be local
            const auto& tuple = procConnections[connecti];
            // const label nbrProci = tuple.first();
            const label nbrFacei = tuple.second();

            inputMeshFaces_[connecti] = nbrFacei;
        }

        return;
    }

    const label nProcs = UPstream::nProcs(comm());
    const label myRank = UPstream::myProcNo(comm());

    const globalIndex globalFaceNum(areaMesh.mesh().nFaces());

    // Boundary inside faces in polyMesh face ids
    const labelList insideFaces
    (
        UIndirectList<label>
        (
            areaMesh.faceLabels(),
            areaMesh.patch().boundaryFaces()
        )
    );


    // Slightly circuitous, but allows maximum reuse of mapDistributeBase

    // 1. Construct a connectivity map using global face numbers

    labelListList connectivity(areaMesh.nBoundaryEdges());
    List<Map<label>> compactMap(nProcs, Map<label>(0));

    // All local mesh faces used
    labelHashSet localUsed(insideFaces);

    forAll(connectivity, connecti)
    {
        labelList& edgeFaces = connectivity[connecti];
        edgeFaces.resize(2);

        // Owner is the boundary inside face (our side)
        // Neighbour is the boundary outside face

        // Connected neighbour
        const auto& tuple = procConnections[connecti];
        const label nbrProci = tuple.first();
        const label nbrFacei = tuple.second();

        if (myRank == nbrProci)
        {
            // Processor-local connectivity
            localUsed.insert(nbrFacei);
        }

        // Global addressing for the connectivity
        edgeFaces[0] = globalFaceNum.toGlobal(insideFaces[connecti]);
        edgeFaces[1] = globalFaceNum.toGlobal(nbrProci, nbrFacei);
    }

    // Create and replace mapping
    static_cast<mapDistributeBase&>(*this) = mapDistributeBase
    (
        globalFaceNum,
        connectivity,
        compactMap,
        UPstream::msgType(),
        comm()
    );

    // List of local mesh faces referenced.
    // Includes inside and locally connected outside faces

    inputMeshFaces_ = localUsed.sortedToc();

    boundaryToCompact_.clear();
    boundaryToCompact_.resize(connectivity.size());

    // After creating the map, connectivity is localized *and*
    // uses compact numbering!

    // Extract the neighbour connection (compact numbering)
    forAll(connectivity, connecti)
    {
        const labelList& edgeFaces = connectivity[connecti];
        // const label face0 = edgeFaces[0];
        const label face1 = edgeFaces[1];

        boundaryToCompact_[connecti] = face1;
    }
}


// ************************************************************************* //
