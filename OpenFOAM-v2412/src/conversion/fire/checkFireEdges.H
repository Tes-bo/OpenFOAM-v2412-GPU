/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

Description
    Checks the mesh for edge connectivity as expected by the AVL/FIRE
    on-the-fly calculations.
    Errors flagged here are not necessarily topological errors at all.

SourceFiles
    checkFireEdges.C

\*---------------------------------------------------------------------------*/

#ifndef checkFireEdges_H
#define checkFireEdges_H

#include "faceList.H"
#include "labelList.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;

//- check edge connectivity
label checkFireEdges
(
    const faceList& faces,
    const labelListList& pointFaces,
    const UList<point>& points = UList<point>::null()
);


//- check edge connectivity with pointFaces mapping calculated automatically
label checkFireEdges
(
    const faceList&,
    const UList<point>& points = UList<point>::null()
);


//- check edge connectivity
label checkFireEdges(const polyMesh& mesh);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
