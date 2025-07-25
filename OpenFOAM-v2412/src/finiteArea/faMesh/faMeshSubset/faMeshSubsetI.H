/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

inline const Foam::faMesh& Foam::faMeshSubset::baseMesh() const noexcept
{
    return baseMesh_;
}


inline const Foam::faMesh& Foam::faMeshSubset::mesh() const noexcept
{
    return subMeshPtr_ ? *subMeshPtr_ : baseMesh_;
}


inline bool Foam::faMeshSubset::hasSubMesh() const noexcept
{
    return bool(subMeshPtr_);
}


inline const Foam::faMesh& Foam::faMeshSubset::subMesh() const
{
    checkHasSubMesh();

    return *subMeshPtr_;
}


inline Foam::faMesh& Foam::faMeshSubset::subMesh()
{
    checkHasSubMesh();

    return *subMeshPtr_;
}


inline const Foam::labelList& Foam::faMeshSubset::pointMap() const
{
    checkHasSubMesh();

    return pointMap_;
}


inline const Foam::labelList& Foam::faMeshSubset::faceMap() const
{
    checkHasSubMesh();

    return faceMap_;
}


inline const Foam::labelList& Foam::faMeshSubset::edgeFlipMap() const
{
    if (!edgeFlipMapPtr_)
    {
        // calcEdgeFlipMap();
    }

    return *edgeFlipMapPtr_;
}


inline const Foam::labelList& Foam::faMeshSubset::cellMap() const
{
    checkHasSubMesh();

    return cellMap_;
}


inline const Foam::labelList& Foam::faMeshSubset::patchMap() const
{
    checkHasSubMesh();

    return patchMap_;
}


// ************************************************************************* //
