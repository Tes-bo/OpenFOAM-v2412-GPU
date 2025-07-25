/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "ensightFaMesh.H"
#include "ensightGeoFile.H"
#include "faMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightFaMesh::clear()
{
    areaPart_.clear();
}


void Foam::ensightFaMesh::renumber()
{
    label partNo = 0;

    areaPart_.index() = partNo++;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFaMesh::ensightFaMesh
(
    const faMesh& mesh
)
:
    mesh_(mesh),
    needsUpdate_(true),
    verbose_(0)
{
    // Lazy?
    if (true)
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::ensightFaMesh::verbose() const noexcept
{
    return verbose_;
}


int Foam::ensightFaMesh::verbose(const int level) noexcept
{
    int old(verbose_);
    verbose_ = level;
    return old;
}


void Foam::ensightFaMesh::correct()
{
    clear();

    // Area meshes (currently only one)
    const label areaId = 0;
    {
        ensightFaces& part = areaPart_;

        part.clear();
        part.identifier() = areaId;
        part.rename("finite-area");

        part.classify
        (
            mesh_.mesh().faces(),
            mesh_.faceLabels()
        );

        // Finalize
        part.reduce();

        if (verbose_)
        {
            Info<< part.info();
        }

        // if (!part.total())
        // {
        //     areaParts_.erase(areaId);
        // }
    }

    renumber();

    needsUpdate_ = false;
}


void Foam::ensightFaMesh::write
(
    ensightGeoFile& os,
    bool parallel
) const
{
    if (UPstream::master())
    {
        os.beginGeometry();
    }

    // Area meshes (currently only one)
    // const label areaId = 0;
    areaPart_.write(os, mesh_.mesh(), parallel);
}


// ************************************************************************* //
