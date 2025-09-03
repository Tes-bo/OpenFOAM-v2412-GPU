/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "normalToFace.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(normalToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, normalToFace, word);
    addToRunTimeSelectionTable(topoSetSource, normalToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, normalToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, normalToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        normalToFace,
        word,
        normal
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        normalToFace,
        istream,
        normal
    );
}


Foam::topoSetSource::addToUsageTable Foam::normalToFace::usage_
(
    normalToFace::typeName,
    "\n    Usage: normalToFace (nx ny nz) <tol>\n\n"
    "    Select faces with normal aligned to unit vector (nx ny nz)\n"
    "    to within tol\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::normalToFace::setNormal()
{
    normal_.normalise();

    if (tol_ < -1 || tol_ > 1)
    {
        FatalErrorInFunction
            << "tolerance not within range -1..1 : " << tol_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalToFace::normalToFace
(
    const polyMesh& mesh,
    const vector& normal,
    const scalar tol
)
:
    topoSetFaceSource(mesh),
    normal_(normal),
    tol_(tol)
{
    setNormal();
}


Foam::normalToFace::normalToFace(const polyMesh& mesh, const dictionary& dict)
:
    topoSetFaceSource(mesh, dict),
    normal_(dict.get<vector>("normal")),
    tol_(dict.get<scalar>("cos"))
{
    setNormal();
}


Foam::normalToFace::normalToFace(const polyMesh& mesh, Istream& is)
:
    topoSetFaceSource(mesh),
    normal_(checkIs(is)),
    tol_(readScalar(checkIs(is)))
{
    setNormal();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::normalToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces according to normal being aligned with "
                << normal_ << " (to within " << tol_ << ") ..." << endl;
        }

        forAll(mesh_.faceAreas(), facei)
        {
            const vector n = normalised(mesh_.faceAreas()[facei]);

            if (mag(1 - (n & normal_)) < tol_)
            {
                set.set(facei);
            }
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces according to normal being aligned with "
                << normal_ << " (to within " << tol_ << ") ..." << endl;
        }

        DynamicList<label> toBeRemoved(set.size()/10);

        for (const label facei : static_cast<const labelHashSet&>(set))
        {
            const vector n = normalised(mesh_.faceAreas()[facei]);

            if (mag(1 - (n & normal_)) < tol_)
            {
                toBeRemoved.append(facei);
            }
        }

        set.unset(toBeRemoved);
    }
}


// ************************************************************************* //
