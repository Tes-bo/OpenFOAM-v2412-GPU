/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label
Foam::ensightFaces::add(const elemType etype, label id, bool flip)
{
    // Linear addressing location
    const label index = offsets_[etype] + sizes_[etype]++;

    addressing()[index] = id;

    if (flipMap_.size())
    {
        flipMap_[index] = flip;
    }

    return index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const char* Foam::ensightFaces::key(const elemType etype) noexcept
{
    return elemNames[etype];
}


inline const Foam::FixedList<Foam::label,3>& Foam::ensightFaces::totals() const
{
    return sizes_;
}


inline Foam::label Foam::ensightFaces::total(const elemType etype) const
{
    return sizes_[etype];
}


inline Foam::label Foam::ensightFaces::size(const elemType etype) const
{
    return (offsets_[etype+1] - offsets_[etype]);
}


inline Foam::labelRange Foam::ensightFaces::range(const elemType etype) const
{
    return labelRange(offsets_[etype], offsets_[etype+1] - offsets_[etype]);
}


inline const Foam::labelList& Foam::ensightFaces::faceIds() const noexcept
{
    return addressing();
}


inline const Foam::labelUList
Foam::ensightFaces::faceIds(const elemType etype) const
{
    return addressing().slice(range(etype));
}


inline const Foam::boolList& Foam::ensightFaces::flipMap() const
{
    return flipMap_;
}


inline bool Foam::ensightFaces::usesFlipMap() const
{
    return (!flipMap_.empty() && flipMap_.size() == size());
}


inline const Foam::labelList&
Foam::ensightFaces::faceOrder() const noexcept
{
    return faceOrder_;
}


inline const Foam::labelUList
Foam::ensightFaces::faceOrder(const elemType etype) const
{
    return faceOrder_.slice(range(etype));
}


inline void Foam::ensightFaces::incrFaceIds(const label off)
{
    incrAddressing(off);
}


inline void Foam::ensightFaces::decrFaceIds(const label off)
{
    decrAddressing(off);
}


// ************************************************************************* //
