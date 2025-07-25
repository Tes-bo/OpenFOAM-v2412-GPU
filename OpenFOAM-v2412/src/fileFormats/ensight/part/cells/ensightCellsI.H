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

inline Foam::label Foam::ensightCells::add(const elemType etype, label id)
{
    // Linear addressing location
    const label index = offsets_[etype] + sizes_[etype]++;

    addressing()[index] = id;

    return index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::ensightCells::manifold() const noexcept
{
    return manifold_;
}


inline const char* Foam::ensightCells::key(const elemType etype) noexcept
{
    return elemNames[etype];
}


inline const Foam::FixedList<Foam::label,5>& Foam::ensightCells::totals() const
{
    return sizes_;
}


inline Foam::label Foam::ensightCells::total(const elemType etype) const
{
    return sizes_[etype];
}


inline Foam::label Foam::ensightCells::size(const elemType etype) const
{
    return (offsets_[etype+1] - offsets_[etype]);
}


inline Foam::labelRange Foam::ensightCells::range(const elemType etype) const
{
    return labelRange(offsets_[etype], offsets_[etype+1] - offsets_[etype]);
}


inline const Foam::labelList& Foam::ensightCells::cellIds() const
{
    return addressing();
}


inline const Foam::labelUList
Foam::ensightCells::cellIds(const elemType etype) const
{
    return addressing().slice(range(etype));
}


inline void Foam::ensightCells::incrCellIds(const label off)
{
    incrAddressing(off);
}


inline void Foam::ensightCells::decrCellIds(const label off)
{
    decrAddressing(off);
}


// ************************************************************************* //
