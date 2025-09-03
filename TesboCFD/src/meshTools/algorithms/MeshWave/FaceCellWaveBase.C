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

#include "FaceCellWave.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FaceCellWaveBase, 0);
}


const Foam::scalar Foam::FaceCellWaveBase::geomTol_ = 1e-6;

Foam::scalar Foam::FaceCellWaveBase::propagationTol_ = 0.01;

int Foam::FaceCellWaveBase::dummyTrackData_ = 12345;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FaceCellWaveBase::FaceCellWaveBase
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    pBufs_(UPstream::commsTypes::nonBlocking),
    changedFace_(mesh_.nFaces()),
    changedCell_(mesh_.nCells()),
    changedFaces_(mesh_.nFaces()),
    changedCells_(mesh_.nCells()),
    nUnvisitedFaces_(mesh_.nFaces()),
    nUnvisitedCells_(mesh_.nCells())
{
    // Don't clear storage on persistent buffer
    pBufs_.allowClearRecv(false);
}


// ************************************************************************* //
