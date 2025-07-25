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

#include "foamVtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::surfaceWriter::surfaceWriter
(
    const pointField& points,
    const faceList& faces,
    const vtk::outputOptions opts
)
:
    vtk::polyWriter(opts),

    points_(std::cref<pointField>(points)),
    faces_(std::cref<faceList>(faces)),
    instant_()
{}


Foam::vtk::surfaceWriter::surfaceWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& file,
    bool parallel
)
:
    surfaceWriter(points, faces)
{
    open(file, parallel);
}


Foam::vtk::surfaceWriter::surfaceWriter
(
    const pointField& points,
    const faceList& faces,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    surfaceWriter(points, faces, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::surfaceWriter::setTime(const instant& inst)
{
    instant_ = inst;
}


bool Foam::vtk::surfaceWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    if (!instant_.name().empty())
    {
        return vtk::fileWriter::beginFile
        (
            "time='" + instant_.name() + "'"
        );
    }

    // Provide default title
    return vtk::fileWriter::beginFile("surface");
}


bool Foam::vtk::surfaceWriter::writeGeometry()
{
    return writePolyGeometry(points_.get(), faces_.get());
}


void Foam::vtk::surfaceWriter::writeTimeValue()
{
    if (!instant_.name().empty())
    {
        vtk::fileWriter::writeTimeValue(instant_.value());
    }
}


void Foam::vtk::surfaceWriter::piece
(
    const pointField& points,
    const faceList& faces
)
{
    endPiece();

    points_ = std::cref<pointField>(points);
    faces_ = std::cref<faceList>(faces);
}


bool Foam::vtk::surfaceWriter::writeProcIDs()
{
    if (this->isPointData())
    {
        return vtk::fileWriter::writeProcIDs(nLocalPoints_);
    }
    return vtk::fileWriter::writeProcIDs(nLocalPolys_);
}


// ************************************************************************* //
