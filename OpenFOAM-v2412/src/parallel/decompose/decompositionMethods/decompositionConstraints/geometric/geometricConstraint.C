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

#include "geometricConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "Time.H"
#include "BitOps.H"
#include "faceBoolSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(geometric);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        geometric,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::geometric::geometric
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    sources_(),
    selection_(coeffDict_.subDict("selection")),
    grow_(dict.getOrDefault("grow", false))
{
    // Stored as dictionary, since we do not have the mesh at this stage

    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding " << selection_.size()
            << " geometric constraints for faces" << endl;
    }
}


Foam::decompositionConstraints::geometric::geometric
(
    PtrList<topoSetFaceSource>&& selections
)
:
    decompositionConstraint(dictionary(), typeName),
    sources_(std::move(selections)),
    selection_(),
    grow_(false)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding " << sources_.size()
            << " geometric constraints for faces" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::geometric::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const label nFaces = mesh.nFaces();

    blockedFace.resize(nFaces, true);

    label nchanged = 0;
    if (decompositionConstraint::debug)
    {
        nchanged = BitOps::count(blockedFace, false);
    }

    // Modify via topoSetFaceSource
    faceBoolSet facesToBlock(mesh, std::move(blockedFace));

    for (const topoSetFaceSource& source : sources_)
    {
        // source.verbose(false);
        source.applyToSet(topoSetSource::SUBTRACT, facesToBlock);
    }

    for (const entry& dEntry : selection_)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Ignoring non-dictionary entry "
                << dEntry << endl;
            continue;
        }

        const dictionary& spec = dEntry.dict();

        auto source = topoSetFaceSource::New
        (
            spec.get<word>("source"),
            mesh,
            spec.optionalSubDict("sourceInfo")
        );
        // source->verbose(false);

        source->applyToSet(topoSetSource::SUBTRACT, facesToBlock);
    }


    // Finished with topo changes
    blockedFace.transfer(facesToBlock.addressing());

    if (decompositionConstraint::debug)
    {
        nchanged = BitOps::count(blockedFace, false) - nchanged;
    }
    else
    {
        nchanged = 0;
    }

    // Grow mode.
    // Include the faces of cells for which there are already two
    // or more faces in a constraint.
    if (grow_)
    {
        bitSet moreUnblocking(nFaces, false);

        label nUnblocked = 0;

        for (label celli=0; celli < mesh.nCells(); ++celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            nUnblocked = 0;
            for (const label facei : cFaces)
            {
                if (!blockedFace[facei])
                {
                    ++nUnblocked;
                    if (nUnblocked > 2)
                    {
                        break;
                    }
                }
            }

            if (nUnblocked > 2)
            {
                moreUnblocking.set(cFaces);
            }
        }

        nUnblocked = 0;

        for (label facei : moreUnblocking)
        {
            if (blockedFace[facei])
            {
                blockedFace[facei] = false;
                ++nUnblocked;
            }
        }

        if (decompositionConstraint::debug)
        {
            Info<< type()
                << " : geometric constraint grow added "
                << returnReduce(nUnblocked, sumOp<label>())
                <<" faces" << endl;
        }

        // Include in the total
        nchanged += nUnblocked;
    }

    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : geometric constraint added for "
            << returnReduce(nchanged, sumOp<label>())
            <<" faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


// ************************************************************************* //
