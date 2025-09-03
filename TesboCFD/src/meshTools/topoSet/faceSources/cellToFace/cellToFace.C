/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "cellToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, word);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, cellToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, cellToFace, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cellToFace::usage_
(
    cellToFace::typeName,
    "\n    Usage: cellToFace <cellSet> all|both|outside\n\n"
    "    Select -all : all faces of cells in the cellSet\n"
    "           -both: faces where both neighbours are in the cellSet\n\n"
);


const Foam::Enum
<
    Foam::cellToFace::cellAction
>
Foam::cellToFace::cellActionNames_
({
    { cellAction::ALL, "all" },
    { cellAction::BOTH, "both" },
    { cellAction::OUTSIDE, "outside" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Selector>
void Foam::cellToFace::combineImpl
(
    topoSet& set,
    const bool add,
    const Selector& cellLabels
) const
{
    if (option_ == ALL)
    {
        // Add all faces from cell
        for (const label celli : cellLabels)
        {
            const labelList& cFaces = mesh_.cells()[celli];

            addOrDelete(set, cFaces, add);
        }
    }
    else if (option_ == BOTH)
    {
        // Add all faces whose both neighbours are in set.

        const label nInt = mesh_.nInternalFaces();
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();


        // Check all internal faces
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (cellLabels.found(own[facei]) && cellLabels.found(nei[facei]))
            {
                addOrDelete(set, facei, add);
            }
        }


        // Get coupled cell status
        boolList neiInSet(mesh_.nBoundaryFaces(), false);

        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    neiInSet[facei-nInt] = cellLabels.found(own[facei]);
                    ++facei;
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiInSet);


        // Check all boundary faces
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    if (cellLabels.found(own[facei]) && neiInSet[facei-nInt])
                    {
                        addOrDelete(set, facei, add);
                    }
                    ++facei;
                }
            }
        }
    }
    else if (option_ == OUTSIDE)
    {
        // Add all faces where only one neighbour is in set.

        const label nInt = mesh_.nInternalFaces();
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();


        // Check all internal faces
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (cellLabels.found(own[facei]) != cellLabels.found(nei[facei]))
            {
                addOrDelete(set, facei, add);
            }
        }


        // Get coupled cell status
        boolList neiInSet(mesh_.nBoundaryFaces(), false);

        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    neiInSet[facei-nInt] = cellLabels.found(own[facei]);
                    ++facei;
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiInSet);


        // Check all boundary faces
        for (const polyPatch& pp : patches)
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                if (cellLabels.found(own[facei]) != neiInSet[facei-nInt])
                {
                    addOrDelete(set, facei, add);
                }
                ++facei;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Selected option is not available"
            << ", option: " << cellActionNames_[option_]
            << exit(FatalError);
    }
}


void Foam::cellToFace::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    if (isZone_)
    {
        const labelList& cellLabels = mesh_.cellZones()[setName];

        combineImpl(set, add, cellLabels);
    }
    else
    {
        // Load the set
        if (!exists(mesh_.time().path()/topoSet::localPath(mesh_, setName)))
        {
            SeriousError
                << "Cannot load set " << setName << endl;
        }

        cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
        const labelHashSet& cellLabels = loadedSet;

        combineImpl(set, add, cellLabels);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const word& setName,
    const cellAction option
)
:
    topoSetFaceSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false),
    option_(option)
{}


Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_)),
    option_(cellActionNames_.get("option", dict))
{}


Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false),
    option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces according to cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, true, setName);
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces according to cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, false, setName);
        }
    }
}


// ************************************************************************* //
