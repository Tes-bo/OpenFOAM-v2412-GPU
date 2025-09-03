/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "clipPlaneToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(clipPlaneToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, clipPlaneToFace, word);
    addToRunTimeSelectionTable(topoSetSource, clipPlaneToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, clipPlaneToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, clipPlaneToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        clipPlaneToFace,
        word,
        clipPlane
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        clipPlaneToFace,
        istream,
        clipPlane
    );
}


Foam::topoSetSource::addToUsageTable Foam::clipPlaneToFace::usage_
(
    clipPlaneToFace::typeName,
    "\n    Usage: clipPlaneToFace (px py pz) (nx ny nz)\n\n"
    "    Select faces with centres above the plane\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::clipPlaneToFace::combine(topoSet& set, const bool add) const
{
    // Face centres above the plane

    const tmp<pointField> tctrs(this->transform(mesh_.faceCentres()));
    const pointField& ctrs = tctrs();

    forAll(ctrs, elemi)
    {
        if (((ctrs[elemi] - point_) & normal_) > 0)
        {
            addOrDelete(set, elemi, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clipPlaneToFace::clipPlaneToFace
(
    const polyMesh& mesh,
    const point& basePoint,
    const vector& normal
)
:
    topoSetFaceSource(mesh),
    point_(basePoint),
    normal_(normal)
{}


Foam::clipPlaneToFace::clipPlaneToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    point_(dict.get<vector>("point")),
    normal_(dict.get<vector>("normal"))
{}


Foam::clipPlaneToFace::clipPlaneToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    point_(checkIs(is)),
    normal_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::clipPlaneToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces with centres above plane at "
                << point_ << " with normal " << normal_ << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces with centres above plane at "
                << point_ << " with normal " << normal_ << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
