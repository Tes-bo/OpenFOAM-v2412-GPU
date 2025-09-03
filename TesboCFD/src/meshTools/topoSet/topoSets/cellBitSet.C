/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "cellBitSet.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(cellBitSet);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellBitSet::cellBitSet(const polyMesh& mesh, const bool val)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), val)
{}


Foam::cellBitSet::cellBitSet
(
    const polyMesh& mesh,
    const bitSet& bits
)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), bits)
{}


Foam::cellBitSet::cellBitSet
(
    const polyMesh& mesh,
    bitSet&& bits
)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), std::move(bits))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellBitSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellBitSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::bitSet Foam::cellBitSet::select
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool verbosity
)
{
    // Start with all cells unselected
    cellBitSet result(mesh);

    // Execute all actions
    for (const entry& dEntry : dict)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Ignoring non-dictionary entry "
                << dEntry << endl;
            continue;
        }

        const dictionary& dict = dEntry.dict();

        const auto action = topoSetSource::combineNames.get("action", dict);

        // These ones we do directly
        switch (action)
        {
            case topoSetSource::INVERT :
            {
                result.invert(mesh.nCells());
                continue;  // Handled
                break;
            }

            case topoSetSource::IGNORE :
                continue;  // Nothing to do
                break;

            default:
                break;
        }

        auto source = topoSetCellSource::New
        (
            dict.get<word>("source"),
            mesh,
            dict.optionalSubDict("sourceInfo")
        );
        source->verbose(verbosity);

        switch (action)
        {
            case topoSetSource::NEW :  // ie, "use"
            case topoSetSource::ADD :
            case topoSetSource::SUBTRACT :
            {
                if (topoSetSource::NEW == action)
                {
                    // "use": only use this selection (CLEAR + ADD)
                    // NEW is handled like ADD in applyToSet()
                    result.reset();
                }
                source->applyToSet(action, result);

                break;
            }

            case topoSetSource::SUBSET :
            {
                cellBitSet other(mesh);
                source->applyToSet(topoSetSource::NEW, other);

                result.subset(other);

                break;
            }

            default:
                // Should already have been caught
                WarningInFunction
                    << "Ignoring unhandled action: "
                    << topoSetSource::combineNames[action] << endl;
        }
    }

    bitSet addr(std::move(result.addressing()));

    return addr;
}


// ************************************************************************* //
