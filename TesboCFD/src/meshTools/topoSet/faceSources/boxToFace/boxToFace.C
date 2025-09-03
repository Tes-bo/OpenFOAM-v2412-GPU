/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "boxToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boxToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, boxToFace, word);
    addToRunTimeSelectionTable(topoSetSource, boxToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, boxToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, boxToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        boxToFace,
        word,
        box
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        boxToFace,
        istream,
        box
    );
}


Foam::topoSetSource::addToUsageTable Foam::boxToFace::usage_
(
    boxToFace::typeName,
    "\n    Usage: boxToFace ((minx miny minz) (maxx maxy maxz))\n\n"
    "    Select all face with faceCentre within bounding box\n\n"
);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read min/max or min/span
static void readBoxDim(const dictionary& dict, treeBoundBox& bb)
{
    IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ;

    dict.readEntry("min", bb.min(), keyType::LITERAL, readOpt);

    if (dict.readIfPresent("span", bb.max(), keyType::LITERAL))
    {
        bb.max() += bb.min();
        readOpt = IOobjectOption::READ_IF_PRESENT;
    }

    dict.readEntry("max", bb.max(), keyType::LITERAL, readOpt);
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boxToFace::combine(topoSet& set, const bool add) const
{
    const tmp<pointField> tctrs(this->transform(mesh_.faceCentres()));
    const pointField& ctrs = tctrs();

    forAll(ctrs, elemi)
    {
        for (const auto& bb : bbs_)
        {
            if (bb.contains(ctrs[elemi]))
            {
                addOrDelete(set, elemi, add);
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    const treeBoundBoxList& bbs
)
:
    topoSetFaceSource(mesh),
    bbs_(bbs)
{}


Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    treeBoundBoxList&& bbs
)
:
    topoSetFaceSource(mesh),
    bbs_(std::move(bbs))
{}


Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    bbs_()
{
    // Accept 'boxes', 'box' or 'min/max'
    if (!dict.readIfPresent("boxes", bbs_))
    {
        bbs_.resize(1);
        if (!dict.readIfPresent("box", bbs_.front()))
        {
            readBoxDim(dict, bbs_.front());
        }
    }
}


Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    bbs_(one{}, treeBoundBox(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boxToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces with centre within boxes "
                << bbs_ << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces with centre within boxes "
                << bbs_ << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
