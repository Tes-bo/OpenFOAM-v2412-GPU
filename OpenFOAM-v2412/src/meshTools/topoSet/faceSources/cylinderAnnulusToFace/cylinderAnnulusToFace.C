/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "cylinderAnnulusToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderAnnulusToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderAnnulusToFace, word);
    addToRunTimeSelectionTable(topoSetSource, cylinderAnnulusToFace, istream);
    addToRunTimeSelectionTable
    (
        topoSetFaceSource,
        cylinderAnnulusToFace,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetFaceSource,
        cylinderAnnulusToFace,
        istream
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        cylinderAnnulusToFace,
        word,
        cylinderAnnulus
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        cylinderAnnulusToFace,
        istream,
        cylinderAnnulus
    );
}


Foam::topoSetSource::addToUsageTable Foam::cylinderAnnulusToFace::usage_
(
    cylinderAnnulusToFace::typeName,
    "\n    Usage: cylinderAnnulusToFace (p1X p1Y p1Z) (p2X p2Y p2Z)"
    " radius innerRadius\n\n"
    "    Select faces with centres within bounding cylinder annulus\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderAnnulusToFace::cylinderAnnulusToFace
(
    const polyMesh& mesh,
    const point& point1,
    const point& point2,
    const scalar radius,
    const scalar innerRadius
)
:
    cylinderToFace(mesh, point1, point2, radius, innerRadius)
{}


Foam::cylinderAnnulusToFace::cylinderAnnulusToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cylinderToFace(mesh, dict)
{}


Foam::cylinderAnnulusToFace::cylinderAnnulusToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    cylinderToFace(mesh, is, true)  // mandatoryInnerRadius = true
{}


// ************************************************************************* //
