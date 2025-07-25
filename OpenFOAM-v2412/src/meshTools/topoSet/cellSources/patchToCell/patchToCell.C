/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "patchToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, patchToCell, word);
    addToRunTimeSelectionTable(topoSetSource, patchToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, patchToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, patchToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        patchToCell,
        word,
        patch
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        patchToCell,
        istream,
        patch
    );
}


Foam::topoSetSource::addToUsageTable Foam::patchToCell::usage_
(
    patchToCell::typeName,
    "\n    Usage: patchToCell patch\n\n"
    "    Select cells attached to patch. Note:accepts wildcards for patch.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToCell::combine(topoSet& set, const bool add) const
{
    labelHashSet patchIDs = mesh_.boundaryMesh().patchSet
    (
        selectedPatches_,
        true,           // warn if not found
        true            // use patch groups if available
    );

    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (verbose_)
        {
            Info<< "    Found matching patch " << pp.name() << " with "
                << returnReduce(pp.size(), sumOp<label>()) << " faces" << endl;
        }

        for
        (
            label facei = pp.start();
            facei < pp.start() + pp.size();
            ++facei
        )
        {
            addOrDelete(set, mesh_.faceOwner()[facei], add);
        }
    }

    if (patchIDs.empty())
    {
        WarningInFunction
            << "Cannot find any patches matching "
            << flatOutput(selectedPatches_) << nl
            << "Valid names: " << flatOutput(mesh_.boundaryMesh().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    const wordRe& patchName
)
:
    topoSetCellSource(mesh),
    selectedPatches_(one{}, patchName)
{}


Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    selectedPatches_()
{
    // Look for 'patches' and 'patch', but accept 'name' as well
    if (!dict.readIfPresent("patches", selectedPatches_))
    {
        selectedPatches_.resize(1);
        selectedPatches_.front() =
            dict.getCompat<wordRe>("patch", {{"name", 1806}});
    }
}


Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    selectedPatches_(one{}, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells associated with patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells associated with patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
