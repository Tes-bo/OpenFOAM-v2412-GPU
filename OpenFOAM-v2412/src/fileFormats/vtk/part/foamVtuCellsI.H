/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "foamVtuCells.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline enum Foam::vtk::vtuCells::contentType
Foam::vtk::vtuCells::content() const noexcept
{
    return output_;
}


inline bool Foam::vtk::vtuCells::decomposeRequested() const noexcept
{
    return decomposeRequest_;
}


inline bool Foam::vtk::vtuCells::empty() const noexcept
{
    return cellTypes_.empty();
}


inline Foam::label Foam::vtk::vtuCells::size() const noexcept
{
    return cellTypes_.size();
}


inline const Foam::List<uint8_t>&
Foam::vtk::vtuCells::cellTypes() const noexcept
{
    return cellTypes_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::vertLabels() const noexcept
{
    return vertLabels_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::vertOffsets() const noexcept
{
    return vertOffset_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::faceLabels() const noexcept
{
    return faceLabels_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::faceOffsets() const noexcept
{
    return faceOffset_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::addPointCellLabels() const noexcept
{
    return maps_.additionalIds();
}


inline const Foam::labelList&
Foam::vtk::vtuCells::cellMap() const noexcept
{
    return maps_.cellMap();
}


// ************************************************************************* //
