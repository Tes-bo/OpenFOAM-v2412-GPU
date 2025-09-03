/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

inline Foam::word Foam::vtk::vtmWriter::ext()
{
    return vtk::fileExtension[vtk::fileTag::MULTI_BLOCK];
}


inline bool Foam::vtk::vtmWriter::append_vtp(const fileName& file)
{
    return append(file, vtk::fileTag::POLY_DATA);
}


inline bool Foam::vtk::vtmWriter::append_vtp
(
    const word& name,
    const fileName& file
)
{
    return append(name, file, vtk::fileTag::POLY_DATA);
}


inline bool Foam::vtk::vtmWriter::append_vtu(const fileName& file)
{
    return append(file, vtk::fileTag::UNSTRUCTURED_GRID);
}


inline bool Foam::vtk::vtmWriter::append_vtu
(
    const word& name,
    const fileName& file
)
{
    return append(name, file, vtk::fileTag::UNSTRUCTURED_GRID);
}


// ************************************************************************* //
