/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "patchToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, patchToFace, word);
    addToRunTimeSelectionTable(topoSetSource, patchToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, patchToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, patchToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        patchToFace,
        word,
        patch
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        patchToFace,
        istream,
        patch
    );
}


Foam::topoSetSource::addToUsageTable Foam::patchToFace::usage_
(
    patchToFace::typeName,
    "\n    Usage: patchToFace patch\n\n"
    "    Select all faces in the patch. Note:accepts wildcards for patch.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToFace::combine(topoSet& set, const bool add) const
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

        for (const label facei : pp.range())
        {
            addOrDelete(set, facei, add);
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

Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    const wordRe& patchName
)
:
    topoSetFaceSource(mesh),
    selectedPatches_(one{}, patchName)
{}


Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
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


Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    selectedPatches_(one{}, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all faces of patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all faces of patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
