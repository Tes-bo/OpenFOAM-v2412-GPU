/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

inline Foam::label Foam::primitiveMesh::nInternalPoints() const noexcept
{
    return nInternalPoints_;
}


inline Foam::label Foam::primitiveMesh::nPoints() const noexcept
{
    return nPoints_;
}


inline Foam::label Foam::primitiveMesh::nInternal0Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal0Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternal1Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal1Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternalEdges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternalEdges_;
}


inline Foam::label Foam::primitiveMesh::nEdges() const
{
    if (nEdges_ < 0)
    {
        nEdges_ = edges().size();
    }

    return nEdges_;
}


inline Foam::label Foam::primitiveMesh::nInternalFaces() const noexcept
{
    return nInternalFaces_;
}


inline Foam::label Foam::primitiveMesh::nBoundaryFaces() const noexcept
{
    return (nFaces_ - nInternalFaces_);
}


inline Foam::label Foam::primitiveMesh::nFaces() const noexcept
{
    return nFaces_;
}


inline Foam::label Foam::primitiveMesh::nCells() const noexcept
{
    return nCells_;
}


inline bool Foam::primitiveMesh::isInternalFace
(
    const label faceIndex
) const noexcept
{
    return faceIndex < nInternalFaces_;
}


inline bool Foam::primitiveMesh::hasCellShapes() const noexcept
{
    return bool(cellShapesPtr_);
}


inline bool Foam::primitiveMesh::hasEdges() const noexcept
{
    return bool(edgesPtr_);
}


inline bool Foam::primitiveMesh::hasCellCells() const noexcept
{
    return bool(ccPtr_);
}


inline bool Foam::primitiveMesh::hasEdgeCells() const noexcept
{
    return bool(ecPtr_);
}


inline bool Foam::primitiveMesh::hasPointCells() const noexcept
{
    return bool(pcPtr_);
}


inline bool Foam::primitiveMesh::hasCells() const noexcept
{
    return bool(cfPtr_);
}


inline bool Foam::primitiveMesh::hasEdgeFaces() const noexcept
{
    return bool(efPtr_);
}


inline bool Foam::primitiveMesh::hasPointFaces() const noexcept
{
    return bool(pfPtr_);
}


inline bool Foam::primitiveMesh::hasCellEdges() const noexcept
{
    return bool(cePtr_);
}


inline bool Foam::primitiveMesh::hasFaceEdges() const noexcept
{
    return bool(fePtr_);
}


inline bool Foam::primitiveMesh::hasPointEdges() const noexcept
{
    return bool(pePtr_);
}


inline bool Foam::primitiveMesh::hasPointPoints() const noexcept
{
    return bool(ppPtr_);
}


inline bool Foam::primitiveMesh::hasCellPoints() const noexcept
{
    return bool(cpPtr_);
}


inline bool Foam::primitiveMesh::hasCellCentres() const noexcept
{
    return bool(cellCentresPtr_);
}


inline bool Foam::primitiveMesh::hasCellVolumes() const noexcept
{
    return bool(cellVolumesPtr_);
}


inline bool Foam::primitiveMesh::hasFaceCentres() const noexcept
{
    return bool(faceCentresPtr_);
}


inline bool Foam::primitiveMesh::hasFaceAreas() const noexcept
{
    return bool(faceAreasPtr_);
}


// ************************************************************************* //
