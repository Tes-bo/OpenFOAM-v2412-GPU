/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "foamVtuSizing.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::vtk::vtuSizing::decompose() const noexcept
{
    return decompose_;
}


inline bool Foam::vtk::vtuSizing::manifold() const noexcept
{
    return manifold_;
}


inline Foam::vtk::vtuSizing::selectionModeType
Foam::vtk::vtuSizing::selectionMode() const noexcept
{
    return selectionMode_;
}


inline Foam::label Foam::vtk::vtuSizing::nCells() const noexcept
{
    return nCells_;
}


inline Foam::label Foam::vtk::vtuSizing::nPoints() const noexcept
{
    return nPoints_;
}


inline Foam::label Foam::vtk::vtuSizing::nVertLabels() const noexcept
{
    return nVertLabels_;
}


inline Foam::label Foam::vtk::vtuSizing::nFaceLabels() const noexcept
{
    return nFaceLabels_;
}


inline Foam::label Foam::vtk::vtuSizing::nCellsPoly() const noexcept
{
    return nCellsPoly_;
}


inline Foam::label Foam::vtk::vtuSizing::nVertPoly() const noexcept
{
    return nVertPoly_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddCells() const noexcept
{
    return nAddCells_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddPoints() const noexcept
{
    return nAddPoints_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddVerts() const noexcept
{
    return nAddVerts_;
}


inline Foam::label Foam::vtk::vtuSizing::nFieldCells() const noexcept
{
    return (nCells_ + nAddCells_);
}


inline Foam::label Foam::vtk::vtuSizing::nFieldPoints() const noexcept
{
    return (nPoints_ + nAddPoints_);
}


inline Foam::label Foam::vtk::vtuSizing::sizeLegacy() const
{
    return sizeOf(contentType::LEGACY, slotType::CELLS);
}


inline Foam::label Foam::vtk::vtuSizing::sizeLegacy
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::LEGACY, slot);
}


inline Foam::label Foam::vtk::vtuSizing::sizeXml
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::XML, slot);
}


inline Foam::label Foam::vtk::vtuSizing::sizeInternal1
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::INTERNAL1, slot);
}


inline Foam::label Foam::vtk::vtuSizing::sizeInternal2
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::INTERNAL2, slot);
}


// ************************************************************************* //
