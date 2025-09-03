/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "PointEdgeWave.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PointEdgeWaveBase, 0);
}


Foam::scalar Foam::PointEdgeWaveBase::propagationTol_ = 0.01;

int Foam::PointEdgeWaveBase::dummyTrackData_ = 12345;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PointEdgeWaveBase::PointEdgeWaveBase
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    pBufs_(UPstream::commsTypes::nonBlocking),
    changedPoint_(mesh_.nPoints()),
    changedEdge_(mesh_.nEdges()),
    changedPoints_(mesh_.nPoints()),
    changedEdges_(mesh_.nEdges()),
    nUnvisitedPoints_(mesh_.nPoints()),
    nUnvisitedEdges_(mesh_.nEdges())
{
    // Don't clear storage on persistent buffer
    pBufs_.allowClearRecv(false);
}


// ************************************************************************* //
