/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField&
Foam::blockDescriptor::vertices() const noexcept
{
    return vertices_;
}


inline const Foam::blockFaceList&
Foam::blockDescriptor::blockFaces() const noexcept
{
    return blockFaces_;
}


inline const Foam::cellShape&
Foam::blockDescriptor::blockShape() const noexcept
{
    return blockShape_;
}


inline const Foam::labelVector&
Foam::blockDescriptor::density() const noexcept
{
    return ijkMesh::sizes();
}


inline const Foam::List<Foam::gradingDescriptors>&
Foam::blockDescriptor::grading() const noexcept
{
    return expand_;
}


inline const Foam::word& Foam::blockDescriptor::zoneName() const noexcept
{
    return zoneName_;
}


inline const Foam::FixedList<Foam::label, 6>&
Foam::blockDescriptor::curvedFaces() const noexcept
{
    return curvedFaces_;
}


inline Foam::label Foam::blockDescriptor::nCurvedFaces() const noexcept
{
    return nCurvedFaces_;
}


inline const Foam::point& Foam::blockDescriptor::blockPoint(const label i) const
{
    return vertices_[blockShape_[i]];
}


inline Foam::label Foam::blockDescriptor::facePointLabel
(
    const direction facei,
    const label i,
    const label j
) const
{
    if (facei == 0 || facei == 1)
    {
        // x-min, x-max
        return
        (
            i
          + j*(sizes().y() + 1)
        );
    }
    else if (facei == 2 || facei == 3)
    {
        // y-min, y-max
        return
        (
            i
          + j*(sizes().x() + 1)
        );
    }
    else
    {
        // z-min, z-max
        return
        (
            i
          + j*(sizes().x() + 1)
        );
    }
}


inline bool Foam::blockDescriptor::vertex
(
    const label i, const label j, const label k
) const
{
    bool iEnd = (i == 0 || i == sizes().x());
    bool jEnd = (j == 0 || j == sizes().y());
    bool kEnd = (k == 0 || k == sizes().z());

    return (iEnd && jEnd && kEnd);
}


inline bool Foam::blockDescriptor::edge
(
    const label i, const label j, const label k
) const
{
    bool iEnd = (i == 0 || i == sizes().x());
    bool jEnd = (j == 0 || j == sizes().y());
    bool kEnd = (k == 0 || k == sizes().z());

    return (iEnd && jEnd) || (iEnd && kEnd) || (jEnd && kEnd);
}


inline bool Foam::blockDescriptor::flatFaceOrEdge
(
    const label i, const label j, const label k
) const
{
    if (i == 0 && curvedFaces_[0] < 0) return true;
    if (i == sizes().x() && curvedFaces_[1] < 0) return true;

    if (j == 0 && curvedFaces_[2] < 0) return true;
    if (j == sizes().y() && curvedFaces_[3] < 0) return true;

    if (k == 0 && curvedFaces_[4] < 0) return true;
    if (k == sizes().z() && curvedFaces_[5] < 0) return true;

    return this->edge(i, j, k);
}


// ************************************************************************* //
