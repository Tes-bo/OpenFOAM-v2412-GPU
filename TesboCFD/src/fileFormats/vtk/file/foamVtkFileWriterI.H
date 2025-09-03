/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline std::ofstream& Foam::vtk::fileWriter::os() noexcept
{
    return os_;
}


inline Foam::vtk::formatter& Foam::vtk::fileWriter::format()
{
    return *format_;
}


inline bool Foam::vtk::fileWriter::isState(outputState test) const noexcept
{
    return (test == state_);
}


inline bool Foam::vtk::fileWriter::notState(outputState test) const noexcept
{
    return (test != state_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::vtk::fileTag Foam::vtk::fileWriter::contentType() const noexcept
{
    return contentType_;
}


inline Foam::vtk::outputOptions Foam::vtk::fileWriter::opts() const noexcept
{
    return opts_;
}


inline Foam::word Foam::vtk::fileWriter::ext() const
{
    return opts_.ext(contentType_);
}


inline bool Foam::vtk::fileWriter::legacy() const noexcept
{
    return opts_.legacy();
}


inline bool Foam::vtk::fileWriter::parallel() const noexcept
{
    return parallel_;
}


inline const Foam::word& Foam::vtk::fileWriter::state() const
{
    return stateNames[state_];
}


inline const Foam::fileName& Foam::vtk::fileWriter::output() const noexcept
{
    return outputFile_;
}


inline Foam::label Foam::vtk::fileWriter::nCellData() const noexcept
{
    return nCellData_;
}


inline Foam::label Foam::vtk::fileWriter::nPointData() const noexcept
{
    return nPointData_;
}


inline bool Foam::vtk::fileWriter::isCellData() const noexcept
{
    return (outputState::CELL_DATA == state_);
}


inline bool Foam::vtk::fileWriter::isPointData() const noexcept
{
    return (outputState::POINT_DATA == state_);
}


// ************************************************************************* //
