/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "foamVtkLineWriter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::lineWriter::lineWriter
(
    const pointField& points,
    const edgeList& edges,
    const vtk::outputOptions opts
)
:
    vtk::polyWriter(opts),

    points_(std::cref<pointField>(points)),
    edges_(std::cref<edgeList>(edges)),
    instant_()
{}


Foam::vtk::lineWriter::lineWriter
(
    const pointField& points,
    const edgeList& edges,
    const fileName& file,
    bool parallel
)
:
    lineWriter(points, edges)
{
    open(file, parallel);
}


Foam::vtk::lineWriter::lineWriter
(
    const pointField& points,
    const edgeList& edges,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    lineWriter(points, edges, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::lineWriter::setTime(const instant& inst)
{
    instant_ = inst;
}


bool Foam::vtk::lineWriter::beginFile(std::string title)
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
    return vtk::fileWriter::beginFile("edges");
}


bool Foam::vtk::lineWriter::writeGeometry()
{
    return writeLineGeometry(points_.get(), edges_.get());
}


void Foam::vtk::lineWriter::writeTimeValue()
{
    if (!instant_.name().empty())
    {
        vtk::fileWriter::writeTimeValue(instant_.value());
    }
}


void Foam::vtk::lineWriter::piece
(
    const pointField& points,
    const edgeList& edges
)
{
    endPiece();

    points_ = std::cref<pointField>(points);
    edges_ = std::cref<edgeList>(edges);
}


bool Foam::vtk::lineWriter::writeProcIDs()
{
    if (this->isPointData())
    {
        return vtk::fileWriter::writeProcIDs(nLocalPoints_);
    }
    return vtk::fileWriter::writeProcIDs(nLocalLines_);
}


// ************************************************************************* //
