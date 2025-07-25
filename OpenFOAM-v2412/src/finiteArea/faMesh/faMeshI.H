/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::faBoundaryMesh& Foam::faMesh::boundary() const noexcept
{
    return boundary_;
}


inline Foam::label Foam::faMesh::nPoints() const noexcept
{
    return nPoints_;
}


inline Foam::label Foam::faMesh::nEdges() const noexcept
{
    return nEdges_;
}


inline Foam::label Foam::faMesh::nInternalEdges() const noexcept
{
    return nInternalEdges_;
}


inline Foam::label Foam::faMesh::nBoundaryEdges() const noexcept
{
    return (nEdges_ - nInternalEdges_);
}


inline Foam::label Foam::faMesh::nFaces() const noexcept
{
    return nFaces_;
}


inline const Foam::pointField& Foam::faMesh::points() const
{
    return patch().localPoints();
}


inline const Foam::edgeList& Foam::faMesh::edges() const noexcept
{
    return edges_;
}


inline const Foam::edgeList::subList Foam::faMesh::internalEdges() const
{
    return edgeList::subList(edges_, nInternalEdges_);
}


inline const Foam::faceList& Foam::faMesh::faces() const
{
    return patch().localFaces();
}


inline const Foam::labelList& Foam::faMesh::edgeOwner() const noexcept
{
    return edgeOwner_;
}


inline const Foam::labelList& Foam::faMesh::edgeNeighbour() const noexcept
{
    return edgeNeighbour_;
}


inline const Foam::labelList& Foam::faMesh::faceLabels() const noexcept
{
    return faceLabels_;
}


inline const Foam::uindirectPrimitivePatch& Foam::faMesh::patch() const
{
    if (!patchPtr_)
    {
        initPatch();
    }
    return *patchPtr_;
}


inline Foam::uindirectPrimitivePatch& Foam::faMesh::patch()
{
    if (!patchPtr_)
    {
        initPatch();
    }
    return *patchPtr_;
}


inline Foam::label Foam::faMesh::whichFace(const label meshFacei) const
{
    return faceLabels_.find(meshFacei);
}


inline const Foam::labelList& Foam::faMesh::whichPolyPatches() const
{
    if (!polyPatchFacesPtr_ || !polyPatchIdsPtr_)
    {
        calcWhichPatchFaces();
    }
    return *polyPatchIdsPtr_;
}


inline const Foam::List<Foam::labelPair>& Foam::faMesh::whichPatchFaces() const
{
    if (!polyPatchFacesPtr_ || !polyPatchIdsPtr_)
    {
        calcWhichPatchFaces();
    }
    return *polyPatchFacesPtr_;
}


inline bool Foam::faMesh::hasInternalEdgeLabels() const noexcept
{
    return false;
}


inline const Foam::List<Foam::labelPair>&
Foam::faMesh::boundaryConnections() const
{
    if (!bndConnectPtr_)
    {
        calcBoundaryConnections();
    }
    return *bndConnectPtr_;
}


// ************************************************************************* //
