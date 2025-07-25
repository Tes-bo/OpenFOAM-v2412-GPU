/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "patchToPoint.H"
#include "pointMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, patchToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, patchToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, patchToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, patchToPoint, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        patchToPoint,
        word,
        patch
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        patchToPoint,
        istream,
        patch
    );
}


Foam::topoSetSource::addToUsageTable Foam::patchToPoint::usage_
(
    patchToPoint::typeName,
    "\n    Usage: patchToPoint patch\n\n"
    "    Select all points in the pointPatch."
    " Note:accepts wildcards for patch.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToPoint::combine(topoSet& set, const bool add) const
{
    const pointMesh& pMesh = pointMesh::New(mesh_, IOobject::READ_IF_PRESENT);
    const pointBoundaryMesh& pbm = pMesh.boundary();

    labelList patchIDs
    (
        pbm.indices
        (
            selectedPatches_,
            true        // useGroups
        )
    );


    for (const label patchi : patchIDs)
    {
        const pointPatch& pp = pbm[patchi];

        if (verbose_)
        {
            Info<< "    Found matching patch " << pp.name() << " with "
                << returnReduce(pp.size(), sumOp<label>()) << " points" << endl;
        }

        for (const label pointi : pp.meshPoints())
        {
            addOrDelete(set, pointi, add);
        }
    }

    if (patchIDs.empty())
    {
        WarningInFunction
            << "Cannot find any patches matching "
            << flatOutput(selectedPatches_) << nl
            //<< "Valid names: " << flatOutput(mesh_.boundaryMesh().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPoint::patchToPoint
(
    const polyMesh& mesh,
    const wordRe& patchName
)
:
    topoSetPointSource(mesh),
    selectedPatches_(one{}, patchName)
{}


Foam::patchToPoint::patchToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh),
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


Foam::patchToPoint::patchToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    selectedPatches_(one{}, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all points of patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all points of patches: "
                << flatOutput(selectedPatches_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
